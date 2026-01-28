using ArgMacros
using CSV
using DataFrames
using BioSequences
using Kmers
using FASTX
using CodecZlib
using ProgressMeter
using Serialization
using CircularArrays

const LOCK = ReentrantLock()

const K = 31

const GENOME_LEN = 367_808

rc(x::Int) = GENOME_LEN - x + 1

function kmerhash(seq)
    hashdict = Dict{UInt, Vector{Int}}()
    for (i,k) in enumerate(FwDNAMers{K}(seq))
        h = fx_hash(k)
        v = haskey(hashdict, h) ? hashdict[h] : Int[]
        hashdict[h] = push!(v, i)
    end
    hashdict
end

struct Mapping
    readrange::UnitRange{Int32}
    strand::Char
    strandrange::UnitRange{Int32}
end

struct Reference
    fkmers::Dict{UInt, Vector{Int}}
    fseqv::CircularArray
    feditsites::Vector{Int}
    rkmers::Dict{UInt, Vector{Int}}
    rseqv::CircularArray
    reditsites::Vector{Int}
end

function extendmatchforward(read, seq, readpos, refpos, editsites)
    readextensionstart = readpos + K - 1
    refextensionstart = refpos + K - 1
    @assert read[readextensionstart] == seq[refextensionstart]
    extension = 0
    while true
        if seq[refextensionstart + extension] ≠ read[readextensionstart + extension]
            (refextensionstart + extension) ∉ editsites && break #not an edit site
            seq[refextensionstart + extension] ≠ DNA_T && break #not edited strand
            read[readextensionstart + extension] ≠ DNA_C && break #only unedited ismatch allowed
        end
        extension += 1
        readpos + K - 1 + extension > length(read) && break
    end
    K + extension - 1
end

function extendmatchbackwards(read, seq, readpos, refpos, editsites)
    readextensionstart = readpos
    refextensionstart = refpos
    @assert read[readextensionstart] == seq[refextensionstart]
    extension = 0
    while true
        if seq[refextensionstart - extension] ≠ read[readextensionstart - extension]
            (refextensionstart - extension) ∉ editsites && break #not an edit site
            seq[refextensionstart - extension] ≠ DNA_T && break #not edited strand
            read[readextensionstart - extension] ≠ DNA_C && break #only unedited ismatch allowed
        end
        extension += 1
        readextensionstart - extension == 0 && break
    end
    K + extension - 1
end

function clockwisedistance(x1, x2, reference)
    x2 > x1 ? x2 - x1 : x2 + length(reference.fseqv) - x1
end

function mapreadtostrand(read, readkmers, ref, startkmerhash, firstpos, strand)
    mappings = Vector{Pair{Mapping, Mapping}}(undef, 0)
    refkmers = strand == '+' ? ref.fkmers : ref.rkmers
    seqv = strand == '+' ? ref.fseqv : ref.rseqv
    editsites = strand == '+' ? ref.feditsites : ref.reditsites
    lastkmer, lastpos = last(collect(readkmers))
    endkmerhash = fx_hash(lastkmer)
    if haskey(refkmers, endkmerhash)
        refstarts = refkmers[startkmerhash]
        refends = refkmers[endkmerhash]
        for s in refstarts, e in refends
            if clockwisedistance(e, s, ref) < 5000
                threeprimematchlength = extendmatchforward(read, seqv, firstpos, s, editsites)
                fiveeprimematchlength = extendmatchbackwards(read, seqv, lastpos, e, editsites)
                if threeprimematchlength + fiveeprimematchlength >= lastpos + K - firstpos - 1
                    push!(mappings, Mapping(firstpos:firstpos+threeprimematchlength-1, strand, s:s+threeprimematchlength-1) =>
                        Mapping(lastpos+K-fiveeprimematchlength:lastpos+K-1, strand, e+K-fiveeprimematchlength:e+K-1))
                end
            end
        end
    end
    mappings
end

function mapread(read::LongDNA{4}, ref::Reference)
    mappings = Vector{Pair{Mapping, Mapping}}(undef, 0)
    readkmers = UnambiguousDNAMers{K}(read)
    isempty(readkmers) && return mappings
    firstkmer, firstpos = first(readkmers)
    startkmerhash = fx_hash(firstkmer)
    if haskey(ref.fkmers, startkmerhash)
        append!(mappings, mapreadtostrand(read, readkmers, ref, startkmerhash, firstpos, '+'))
    end
    if haskey(ref.rkmers, startkmerhash)
        append!(mappings, mapreadtostrand(read, readkmers, ref, startkmerhash, firstpos, '-'))
    end
    mappings
end

function main(ARGS)

    @inlinearguments begin
        @positionalrequired String readsfile
        @positionalrequired String outprefix
        @positionalrequired String references
    end

    #read reference
    id = first(split(basename(ARGS[3]), "."))
    gff = CSV.File("$id.gff"; comment = "#", header = ["accession", "software", "feature", "start", "stop", "score", "strand", "phase", "attributes"]) |> DataFrame
    editsites = filter(x -> occursin("C to U RNA editing", x.attributes), gff)
    reference = FASTA.Reader(open(ARGS[3])) do infile
        fseq = FASTA.sequence(LongDNA{2}, first(infile))
        #construct reference kmer arrays
        fkmers = kmerhash(fseq * fseq[1:30])
        fseqv = CircularArray(collect(fseq))
        feditsites = filter(x -> fseqv[x.start] == DNA_T, editsites)
        rseq = FASTA.sequence(LongDNA{2}, first(infile))
        rkmers = kmerhash(rseq * rseq[1:30])
        rseqv = CircularArray(collect(rseq))
        reditsites = filter(x -> rseqv[rc(x.start)] == DNA_T, editsites)
        Reference(fkmers, fseqv, feditsites.start, rkmers, rseqv, reditsites.start)
    end

    #read reads
    reader = FASTQ.Reader(GzipDecompressorStream(open(readsfile)))

    prog = ProgressUnknown("Reads: ")
    numthreads = Threads.nthreads()
    reads = [FASTQ.Record() for i=1:numthreads]
    circularmappings = Vector{Pair{Mapping, Mapping}}(undef,0)
    while !eof(reader)
        for t in 1:numthreads
            read!(reader, reads[t])
            if eof(reader)
                resize!(reads, t)
                break
            end
        end
        Threads.@threads for r in 1:length(reads)
            mappings = mapread(sequence(LongDNA{4}, reads[r]) , reference)
            if ~isempty(mappings)
                lock(LOCK)
                append!(circularmappings, mappings)
                unlock(LOCK)
            end
            ProgressMeter.next!(prog)
        end
    end
    close(reader)
    ProgressMeter.finish!(prog)

    #write results
    serialize(outprefix * ".mappings3.bin", circularmappings)

end
@main

#ARGS[1] = reads (.fq.gz)
#ARGS[2] = output prefix
#ARGS[3] = reference(s) (.fa)

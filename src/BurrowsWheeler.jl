module BurrowsWheeler
using BioSequences
export bwt_naïve, bwt, sa, rank
    """
        sa(word::BioSequence{DNAAlphabet{4}})::Array{Int64}

    Returns the indices of lexicologically sorted suffix array from DNASequence.
    """
    function sa(word::BioSequence{DNAAlphabet{4}})
        w = copy(word)
        if DNA_Gap in w
            throw(ArgumentError("input must not contain gaps"))
        end
        push!(w, DNA_Gap)
        S = sort!([(i, w[i:end]) for i in 1:(length(w))], by = x -> x[2])
        return(first.(S))
    end
    """
        bwt(word::BioSequence{DNAAlphabet{4}})::Array{DNA}

    Returns BurrowsWheeler transformation using sorted suffix array.
    """
    function bwt(word::BioSequence{DNAAlphabet{4}})
        l = DNA[]
        if DNA_Gap in word
            throw(ArgumentError("input must not contain gaps"))
        end
        for i in sa(word)
            j = i - 1
            if j == 0
                push!(l, DNA_Gap)
            else
                push!(l, word[j])
            end
        end
        return(l)
    end
    """
        bwt_naïve(word::BioSequence{DNAAlphabet{4}})::Array{DNA}

    Returns BurrowsWheeler transformation using naïve implementation with two dimensional array.
    """
    function bwt_naïve(word::BioSequence{DNAAlphabet{4}})
        if DNA_Gap in word
            throw(ArgumentError("input must not contain gaps"))
        end
        w = convert(Vector{DNA}, word)
        push!(w, DNA_Gap)
        l = Vector{DNA}()
        i = length(w)
        while i > 0
            c = pop!(w)
            pushfirst!(w, c)
            i = i - 1
            append!(l, w)
        end
        M = reshape(l, length(w), length(w))
        return(sortslices(M, dims=2)[end, :])
    end

    """
        rank(bw::Array{DNA})::Dict{DNA,Int64}

    Returns ranks for the bwt(W) and counts for occurences of each symbol .
    """
    function rank(bw::Array{DNA})
        counts = Dict{DNA,Int64}()
        ranks = Int64[]
        for c in bw
            if !haskey(counts, c)
                counts[c] = 0
            end
            push!(ranks, counts[c])
            counts[c] = counts[c] + 1
        end
        return (ranks, counts)
    end
end

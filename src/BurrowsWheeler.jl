module BurrowsWheeler
using BioSequences
export bwt_naïve, bwt, sa, rank, bwmcol, revbwt
    """
        sa(word::BioSequence{DNAAlphabet{4}})

    Returns the indices of lexicographically sorted suffix array from DNASequence.
    """
    function sa(word::BioSequence{DNAAlphabet{4}})
        w = copy(word)
        if DNA_Gap in w
            throw(ArgumentError("input must not contain gaps"))
        end
        push!(w, DNA_Gap)
        S = sort!([(i, w[i:end]) for i in 1:length(w)], by = x -> x[2])
        return first.(S)::Vector{Int64}
    end
    """
        bwt(word::BioSequence{DNAAlphabet{4}})

    Returns BurrowsWheeler transformation using sorted suffix array.
    """
    function bwt(word::BioSequence{DNAAlphabet{4}})
        if DNA_Gap in word
            throw(ArgumentError("input must not contain gaps"))
        end
        sa_result = sa(word)
        l = Vector{DNA}(undef, length(sa_result))
        @inbounds for (idx, i) in enumerate(sa_result)
            j = i - 1
            if j == 0
                l[idx] = DNA_Gap
            else
                l[idx] = word[j]
            end
        end
        return l
    end
    """
        bwt_naïve(word::BioSequence{DNAAlphabet{4}})::Array{DNA}

    Returns BurrowsWheeler transformation using naïve implementation with two dimensional array.
    """
    function bwt_naïve(word::BioSequence{DNAAlphabet{4}})
        if DNA_Gap in word
            throw(ArgumentError("input must not contain gaps"))
        end
        w = collect(word)
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
        rank(bw::Vector{DNA})

    Returns ranks for the bwt(W) and counts for occurrences of each symbol.
    """
    function rank(bw::Vector{DNA})
        counts = Dict{DNA,Int64}()
        ranks = Vector{Int64}(undef, length(bw))
        for (idx, c) in enumerate(bw)
            count = get!(counts, c, 0)
            ranks[idx] = count
            counts[c] = count + 1
        end
        return (ranks, counts)
    end

    """
        bwmcol(counts::Dict{DNA,Int64})

    Returns counts for the first column of bwt_naïve().
    """
    function bwmcol(counts::Dict{DNA,Int64})
        col = Dict{DNA, Tuple{Int64, Int64}}()
        sums = 0
        for (c, count) in sort(collect(counts), by=x->x[1])
            col[c] = (sums, sums + count)
            sums = sums + count
        end
        return col
    end

    """
        revbwt(bw::Vector{DNA})

    Returns reverse of Burrows Wheeler transformation bwt() function.
    """
    function revbwt(bw::Vector{DNA})
        ranks, counts = rank(bw)
        col = bwmcol(counts)
        rowi = 1
        result = Vector{DNA}(undef, length(bw))
        result[1] = DNA_Gap
        result_len = 1
        @inbounds while bw[rowi] != DNA_Gap
            c = bw[rowi]
            result_len += 1
            result[result_len] = c
            rowi = col[c][1] + ranks[rowi] + 1
        end
        reverse!(result, 1, result_len)
        return LongSequence{DNAAlphabet{4}}(result[1:result_len-1])
    end
end

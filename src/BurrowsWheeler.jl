module BurrowsWheeler
using BioSequences
export bwt_naïve, bwt, sa, rank, bwmcol, revbwt

    _alphabet_type(::Type{DNA}) = DNAAlphabet{4}
    _alphabet_type(::Type{RNA}) = RNAAlphabet{4}

    """
        sa(word::BioSequence{A}) where {A <: Alphabet}

    Returns the indices of lexicographically sorted suffix array.
    """
    function sa(word::BioSequence{A}) where {A <: Alphabet}
        w = copy(word)
        gap_sym = gap(eltype(w))
        if gap_sym in w
            throw(ArgumentError("input must not contain gaps"))
        end
        push!(w, gap_sym)
        S = sort!([(i, w[i:end]) for i in 1:length(w)], by = x -> x[2])
        return first.(S)::Vector{Int64}
    end
    """
        bwt(word::BioSequence{A}) where {A <: Alphabet}

    Returns BurrowsWheeler transformation using sorted suffix array.
    """
    function bwt(word::BioSequence{A}) where {A <: Alphabet}
        S = eltype(word)
        gap_sym = gap(S)
        if gap_sym in word
            throw(ArgumentError("input must not contain gaps"))
        end
        sa_result = sa(word)
        l = Vector{S}(undef, length(sa_result))
        @inbounds for (idx, i) in enumerate(sa_result)
            j = i - 1
            if j == 0
                l[idx] = gap_sym
            else
                l[idx] = word[j]
            end
        end
        return l
    end
    """
        bwt_naïve(word::BioSequence{A}) where {A <: Alphabet}

    Returns BurrowsWheeler transformation using naïve implementation with two dimensional array.
    """
    function bwt_naïve(word::BioSequence{A}) where {A <: Alphabet}
        S = eltype(word)
        gap_sym = gap(S)
        if gap_sym in word
            throw(ArgumentError("input must not contain gaps"))
        end
        w = collect(word)
        push!(w, gap_sym)
        l = Vector{S}()
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
        rank(bw::Vector{S}) where {S}

    Returns ranks for the bwt(W) and counts for occurrences of each symbol.
    """
    function rank(bw::Vector{S}) where {S}
        counts = Dict{S,Int64}()
        ranks = Vector{Int64}(undef, length(bw))
        for (idx, c) in enumerate(bw)
            count = get!(counts, c, 0)
            ranks[idx] = count
            counts[c] = count + 1
        end
        return (ranks, counts)
    end

    """
        bwmcol(counts::Dict{S,Int64}) where {S}

    Returns counts for the first column of bwt_naïve().
    """
    function bwmcol(counts::Dict{S,Int64}) where {S}
        col = Dict{S, Tuple{Int64, Int64}}()
        sums = 0
        for (c, count) in sort(collect(counts), by=x->x[1])
            col[c] = (sums, sums + count)
            sums = sums + count
        end
        return col
    end

    """
        revbwt(bw::Vector{S}, alphabet_type::Type{A}) where {S, A <: Alphabet}

    Returns reverse of Burrows Wheeler transformation bwt() function.
    """
    function revbwt(bw::Vector{S}, alphabet_type::Type{A}) where {S, A <: Alphabet}
        ranks, counts = rank(bw)
        col = bwmcol(counts)
        gap_sym = gap(S)
        rowi = 1
        result = Vector{S}(undef, length(bw))
        result[1] = gap_sym
        result_len = 1
        @inbounds while bw[rowi] != gap_sym
            c = bw[rowi]
            result_len += 1
            result[result_len] = c
            rowi = col[c][1] + ranks[rowi] + 1
        end
        reverse!(result, 1, result_len)
        return LongSequence{A}(result[1:result_len-1])
    end

    """
        revbwt(bw::Vector{S}) where {S}

    Returns reverse of Burrows Wheeler transformation bwt() function.
    """
    function revbwt(bw::Vector{S}) where {S}
        ranks, counts = rank(bw)
        col = bwmcol(counts)
        gap_sym = gap(S)
        rowi = 1
        result = Vector{S}(undef, length(bw))
        result[1] = gap_sym
        result_len = 1
        @inbounds while bw[rowi] != gap_sym
            c = bw[rowi]
            result_len += 1
            result[result_len] = c
            rowi = col[c][1] + ranks[rowi] + 1
        end
        reverse!(result, 1, result_len)
        A = _alphabet_type(S)
        return LongSequence{A}(result[1:result_len-1])
    end
end

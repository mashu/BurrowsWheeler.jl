module Search

using BioSequences
using SuffixArrays
using ..BurrowsWheeler

export FMIndex, build_fmindex, search, count_occurrences, locate

struct FMIndex{S, A <: Alphabet}
    bwt::Vector{S}
    C::Dict{S, Int64}
    rank_arrays::Dict{S, Vector{Int64}}
    sa_sample::Vector{Int64}
    sa_sample_pos::Vector{Int64}
    sample_rate::Int64
end

function _build_rank_arrays(bwt::Vector{S}) where {S}
    n = length(bwt)
    rank_arrays = Dict{S, Vector{Int64}}()
    
    unique_symbols = unique(bwt)
    running_counts = Dict{S, Int64}()
    for sym in unique_symbols
        rank_arrays[sym] = Vector{Int64}(undef, n)
        running_counts[sym] = 0
    end
    
    @inbounds for idx in 1:n
        c = bwt[idx]
        running_counts[c] += 1
        for sym in unique_symbols
            rank_arrays[sym][idx] = running_counts[sym]
        end
    end
    
    return rank_arrays
end

function _build_C(counts::Dict{S, Int64}) where {S}
    C = Dict{S, Int64}()
    sums = 0
    for (c, count) in sort(collect(counts), by=x->x[1])
        C[c] = sums
        sums += count
    end
    return C
end

function _rank(bwt::Vector{S}, c::S, pos::Int64, rank_arrays::Dict{S, Vector{Int64}}) where {S}
    if pos < 1
        return 0
    end
    if pos > length(bwt)
        pos = length(bwt)
    end
    @inbounds return rank_arrays[c][pos]
end

function backward_search(fm::FMIndex{S, A}, pattern::BioSequence{A}) where {S, A <: Alphabet}
    sp = 1
    ep = length(fm.bwt)
    
    @inbounds for i in length(pattern):-1:1
        c = pattern[i]
        
        c_val = get(fm.C, c, -1)
        if c_val == -1
            return (0, 0)
        end
        
        sp = c_val + _rank(fm.bwt, c, sp - 1, fm.rank_arrays) + 1
        ep = c_val + _rank(fm.bwt, c, ep, fm.rank_arrays)
        
        if sp > ep
            return (0, 0)
        end
    end
    
    return (sp, ep)
end

function _sample_suffix_array(sa::Vector{Int64}, sample_rate::Int64)
    sample_pos = Int64[]
    sample_val = Int64[]
    
    @inbounds for (idx, pos) in enumerate(sa)
        if pos == 1 || (pos - 1) % sample_rate == 0
            push!(sample_pos, idx)
            push!(sample_val, pos)
        end
    end
    
    return (sample_pos, sample_val)
end

function _locate_row(fm::FMIndex{S, A}, row::Int64) where {S, A <: Alphabet}
    if row < 1 || row > length(fm.bwt)
        return 0
    end
    
    steps = 0
    current_row = row
    
    while current_row > 0 && current_row <= length(fm.bwt)
        idx = findfirst(x -> x == current_row, fm.sa_sample_pos)
        if idx !== nothing
            @inbounds return fm.sa_sample[idx] + steps
        end
        
        @inbounds c = fm.bwt[current_row]
        rank_arr = get(fm.rank_arrays, c, Int64[])
        if isempty(rank_arr)
            return 0
        end
        @inbounds rank_val = rank_arr[current_row]
        c_val = get(fm.C, c, 0)
        current_row = c_val + rank_val
        steps += 1
        
        if steps > length(fm.bwt)
            return 0
        end
    end
    
    return 0
end

function build_fmindex(sequence::BioSequence{A}; sample_rate::Int64=32) where {A <: Alphabet}
    S = eltype(sequence)
    
    gap_sym = gap(S)
    if gap_sym in sequence
        throw(ArgumentError("input must not contain gaps"))
    end
    
    encoded = BurrowsWheeler._encode_sequence(sequence)
    push!(encoded, 0x00)
    
    sa_0based = suffixsort(encoded, 0)
    sa_result = Vector{Int64}(undef, length(sa_0based))
    @inbounds for (idx, s) in enumerate(sa_0based)
        sa_result[idx] = Int(s) + 1
    end
    
    bwt_result = Vector{S}(undef, length(sa_result))
    @inbounds for (idx, i) in enumerate(sa_result)
        j = i - 1
        if j == 0
            bwt_result[idx] = gap_sym
        else
            bwt_result[idx] = sequence[j]
        end
    end
    
    counts = Dict{S,Int64}()
    ranks = Vector{Int64}(undef, length(bwt_result))
    for (idx, c) in enumerate(bwt_result)
        count = get!(counts, c, 0)
        ranks[idx] = count
        counts[c] = count + 1
    end
    
    C = _build_C(counts)
    rank_arrays = _build_rank_arrays(bwt_result)
    
    sa_sample_pos, sa_sample_val = _sample_suffix_array(sa_result, sample_rate)
    
    return FMIndex{S, A}(
        bwt_result,
        C,
        rank_arrays,
        sa_sample_val,
        sa_sample_pos,
        sample_rate
    )
end

function search(fm::FMIndex{S, A}, pattern::BioSequence{A}) where {S, A <: Alphabet}
    return backward_search(fm, pattern)
end

function count_occurrences(fm::FMIndex{S, A}, pattern::BioSequence{A}) where {S, A <: Alphabet}
    sp, ep = backward_search(fm, pattern)
    if sp < 1 || sp > ep
        return 0
    end
    return ep - sp + 1
end

function locate(fm::FMIndex{S, A}, pattern::BioSequence{A}) where {S, A <: Alphabet}
    sp, ep = backward_search(fm, pattern)
    if sp < 1 || sp > ep
        return Int64[]
    end
    
    positions = Vector{Int64}(undef, ep - sp + 1)
    @inbounds for (idx, row) in enumerate(sp:ep)
        positions[idx] = _locate_row(fm, row)
    end
    
    sort!(positions)
    return positions
end

end


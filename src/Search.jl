module Search

using BioSequences
using SuffixArrays
using ..BurrowsWheeler

export FMIndex, build_fmindex, search, count_occurrences, locate
export ApproxMatch, approximate_search, approximate_locate

"""
    FMIndex{S, A <: Alphabet}

FM-index data structure for efficient pattern matching in genomic sequences.

# Fields
- `bwt`: Burrows-Wheeler Transform of the sequence
- `C`: First column mapping (cumulative counts)
- `rank_arrays`: Rank arrays for each symbol
- `sa_sample`: Sampled suffix array values
- `sa_sample_pos`: Positions of sampled suffix array entries
- `sample_rate`: Sampling rate for suffix array
"""
struct FMIndex{S, A <: Alphabet}
    bwt::Vector{S}
    C::Dict{S, Int64}
    rank_arrays::Dict{S, Vector{Int64}}
    sa_sample::Vector{Int64}
    sa_sample_pos::Vector{Int64}
    sample_rate::Int64
end

"""
    ApproxMatch

Represents an approximate match found during approximate search.

# Fields
- `position`: 1-based position in the original sequence where the match starts
- `edit_distance`: Edit distance (number of substitutions, insertions, deletions)
- `alignment_length`: Length of the alignment (may differ from pattern length due to indels)
"""
struct ApproxMatch
    position::Int64
    edit_distance::Int64
    alignment_length::Int64
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

"""
Backward search algorithm for FM-index pattern matching.

This function implements the backward search algorithm which processes the pattern from right to left,
narrowing the suffix array range [sp, ep] at each step. All suffixes in the final range [sp, ep]
begin with the pattern, ensuring they represent valid substring matches in the original sequence.
"""
function backward_search(fm::FMIndex{S, A}, pattern::BioSequence{A}) where {S, A <: Alphabet}
    sp = 1
    ep = length(fm.bwt)
    
    # Process pattern from right to left to narrow the suffix array range
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

"""
    build_fmindex(sequence::BioSequence{A}; sample_rate::Int64=32) where {A <: Alphabet}

Build an FM-index from a genomic sequence for efficient pattern matching.

# Arguments
- `sequence`: The genomic sequence to index (DNA or RNA)
- `sample_rate`: Sampling rate for suffix array (default: 32). Lower values use more memory but enable faster position reporting.

# Returns
An `FMIndex` structure that can be used for pattern searching.

# Examples
```julia
using BurrowsWheeler
using BioSequences

sequence = dna"ATAATA"
fm = build_fmindex(sequence)
```
"""
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

"""
    search(fm::FMIndex{S, A}, pattern::BioSequence{A}) where {S, A <: Alphabet}

Search for a pattern in the FM-index using backward search.

The backward search algorithm finds all suffixes in the suffix array that begin with the pattern.
This ensures that all positions in the returned range represent valid substring matches in the
original sequence. The algorithm processes the pattern from right to left, narrowing the suffix
array range at each step to guarantee correctness.

# Arguments
- `fm`: The FM-index to search in
- `pattern`: The pattern to search for

# Returns
A tuple `(sp, ep)` where `sp` is the start position and `ep` is the end position in the suffix array.
All suffixes in the range [sp, ep] begin with the pattern, representing valid substring matches.
If the pattern is not found, returns `(0, 0)`.

# Examples
```julia
fm = build_fmindex(dna"ATAATA")
sp, ep = search(fm, dna"ATA")
```
"""
function search(fm::FMIndex{S, A}, pattern::BioSequence{A}) where {S, A <: Alphabet}
    return backward_search(fm, pattern)
end

"""
    count_occurrences(fm::FMIndex{S, A}, pattern::BioSequence{A}) where {S, A <: Alphabet}

Count the number of occurrences of a pattern in the indexed sequence.

# Arguments
- `fm`: The FM-index to search in
- `pattern`: The pattern to count

# Returns
The number of occurrences of the pattern (0 if not found).

# Examples
```julia
fm = build_fmindex(dna"ATAATA")
count = count_occurrences(fm, dna"ATA")  # Returns 2
```
"""
function count_occurrences(fm::FMIndex{S, A}, pattern::BioSequence{A}) where {S, A <: Alphabet}
    sp, ep = backward_search(fm, pattern)
    if sp < 1 || sp > ep
        return 0
    end
    return ep - sp + 1
end

"""
    locate(fm::FMIndex{S, A}, pattern::BioSequence{A}) where {S, A <: Alphabet}

Find all positions where a pattern occurs in the original sequence.

This function uses backward search to find all suffixes beginning with the pattern, then converts
suffix array positions to actual sequence positions. Each returned position represents a valid
substring match: the substring starting at that position in the original sequence matches the
pattern exactly.

# Arguments
- `fm`: The FM-index to search in
- `pattern`: The pattern to locate

# Returns
A sorted vector of 1-based positions where the pattern occurs in the original sequence.
Each position `pos` represents a substring match: `sequence[pos:pos+length(pattern)-1] == pattern`.
Returns an empty vector if the pattern is not found.

# Examples
```julia
fm = build_fmindex(dna"ATAATA")
positions = locate(fm, dna"ATA")  # Returns [1, 4]
# Both positions represent valid matches:
# sequence[1:3] == dna"ATA" and sequence[4:6] == dna"ATA"
```
"""
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

"""
    _backward_search_step(fm::FMIndex{S, A}, c::S, sp::Int64, ep::Int64) where {S, A <: Alphabet}

Perform a single step of backward search for character `c` in the range [sp, ep].

Returns the new range (sp_new, ep_new) after extending with character `c`.
"""
function _backward_search_step(fm::FMIndex{S, A}, c::S, sp::Int64, ep::Int64) where {S, A <: Alphabet}
    c_val = get(fm.C, c, -1)
    if c_val == -1
        return (0, 0)
    end
    
    sp_new = c_val + _rank(fm.bwt, c, sp - 1, fm.rank_arrays) + 1
    ep_new = c_val + _rank(fm.bwt, c, ep, fm.rank_arrays)
    
    if sp_new > ep_new
        return (0, 0)
    end
    
    return (sp_new, ep_new)
end

"""
    _approximate_search_backtrack(fm::FMIndex{S, A}, pattern::BioSequence{A}, 
                                   sp::Int64, ep::Int64, pattern_idx::Int64,
                                   max_edits::Int64, current_edits::Int64,
                                   text_chars_consumed::Int64,
                                   results::Vector{ApproxMatch}, 
                                   visited::Set{Tuple{Int64, Int64, Int64}}) where {S, A <: Alphabet}

Recursive backtracking function for approximate search with edit distance.

This function explores all possible paths in the FM-index that match the pattern within
the maximum edit distance, allowing substitutions, insertions, and deletions.

# Arguments
- `text_chars_consumed`: Number of characters from the text that have been consumed so far.
                         This is used to calculate the actual alignment length.
"""
function _approximate_search_backtrack(fm::FMIndex{S, A}, pattern::BioSequence{A}, 
                                       sp::Int64, ep::Int64, pattern_idx::Int64,
                                       max_edits::Int64, current_edits::Int64,
                                       text_chars_consumed::Int64,
                                       results::Vector{ApproxMatch}, 
                                       visited::Set{Tuple{Int64, Int64, Int64}}) where {S, A <: Alphabet}
    
    # Check if we've already visited this state
    state = (sp, ep, pattern_idx)
    if state in visited
        return
    end
    push!(visited, state)
    
    # Base case: pattern is fully consumed
    if pattern_idx < 1
        # We found a match! Add all positions in the range [sp, ep]
        # alignment_length = text_chars_consumed (actual number of text chars matched)
        if sp >= 1 && sp <= ep
            for row in sp:ep
                pos = _locate_row(fm, row)
                if pos > 0
                    push!(results, ApproxMatch(pos, current_edits, text_chars_consumed))
                end
            end
        end
        return
    end
    
    # If we've exceeded max edits, prune this branch
    if current_edits > max_edits
        return
    end
    
    c = pattern[pattern_idx]
    
    # Try exact match (no edit)
    # Consumes one text character and advances pattern by one
    sp_exact, ep_exact = _backward_search_step(fm, c, sp, ep)
    if sp_exact <= ep_exact
        _approximate_search_backtrack(fm, pattern, sp_exact, ep_exact, pattern_idx - 1,
                                     max_edits, current_edits, text_chars_consumed + 1,
                                     results, visited)
    end
    
    # If we can still make edits, try all possible edits
    if current_edits < max_edits
        # Try substitution: match any character except c
        # Consumes one text character and advances pattern by one
        for sym in keys(fm.C)
            if sym != c
                sp_sub, ep_sub = _backward_search_step(fm, sym, sp, ep)
                if sp_sub <= ep_sub
                    _approximate_search_backtrack(fm, pattern, sp_sub, ep_sub, pattern_idx - 1,
                                                 max_edits, current_edits + 1, text_chars_consumed + 1,
                                                 results, visited)
                end
            end
        end
        
        # Try insertion: consume a character from the text without advancing in pattern
        # This means we stay at pattern_idx but extend the range (text_chars_consumed increases)
        for sym in keys(fm.C)
            sp_ins, ep_ins = _backward_search_step(fm, sym, sp, ep)
            if sp_ins <= ep_ins
                _approximate_search_backtrack(fm, pattern, sp_ins, ep_ins, pattern_idx,
                                             max_edits, current_edits + 1, text_chars_consumed + 1,
                                             results, visited)
            end
        end
        
        # Try deletion: advance in pattern without consuming a character from text
        # This means we decrease pattern_idx but keep the same range (text_chars_consumed stays same)
        _approximate_search_backtrack(fm, pattern, sp, ep, pattern_idx - 1,
                                     max_edits, current_edits + 1, text_chars_consumed,
                                     results, visited)
    end
end

"""
    approximate_search(fm::FMIndex{S, A}, pattern::BioSequence{A}, 
                      max_edits::Int64) where {S, A <: Alphabet}

Find all approximate matches of a pattern in the FM-index within a maximum edit distance.

This function uses a backtracking algorithm to explore all possible alignments allowing
substitutions, insertions, and deletions. The edit distance is the total number of
these operations needed to align the pattern with a substring of the text.

# Arguments
- `fm`: The FM-index to search in
- `pattern`: The pattern to search for
- `max_edits`: Maximum allowed edit distance (default: 2)

# Returns
A vector of `ApproxMatch` objects, each containing:
- `position`: 1-based position in the original sequence
- `edit_distance`: Number of edits (substitutions, insertions, deletions)
- `alignment_length`: Length of the alignment

# Examples
```julia
fm = build_fmindex(dna"ATCGATCG")
matches = approximate_search(fm, dna"ATGG", max_edits=1)
# Finds matches with up to 1 edit (e.g., ATCG at position 1 with substitution)
```
"""
function approximate_search(fm::FMIndex{S, A}, pattern::BioSequence{A}; 
                            max_edits::Int64=2) where {S, A <: Alphabet}
    results = ApproxMatch[]
    visited = Set{Tuple{Int64, Int64, Int64}}()
    
    sp = 1
    ep = length(fm.bwt)
    
    # Start with text_chars_consumed = 0 (no characters consumed yet)
    _approximate_search_backtrack(fm, pattern, sp, ep, length(pattern),
                                 max_edits, 0, 0, results, visited)
    
    # Remove duplicates and sort by position, then by edit distance
    unique_results = Dict{Int64, ApproxMatch}()
    for match in results
        key = match.position
        if !haskey(unique_results, key) || match.edit_distance < unique_results[key].edit_distance
            unique_results[key] = match
        end
    end
    
    sorted_results = sort(collect(values(unique_results)), 
                          by=x -> (x.position, x.edit_distance))
    return sorted_results
end

"""
    approximate_locate(fm::FMIndex{S, A}, pattern::BioSequence{A}, 
                      max_edits::Int64) where {S, A <: Alphabet}

Find the best approximate match positions for a pattern within a maximum edit distance.

This is a convenience function that returns just the positions of the best matches
(smallest edit distance) for each unique position.

# Arguments
- `fm`: The FM-index to search in
- `pattern`: The pattern to search for
- `max_edits`: Maximum allowed edit distance (default: 2)

# Returns
A vector of 1-based positions representing the best approximate matches.
Each position appears at most once, and matches are sorted by position.

# Examples
```julia
fm = build_fmindex(dna"ATCGATCG")
positions = approximate_locate(fm, dna"ATGG", max_edits=1)
```
"""
function approximate_locate(fm::FMIndex{S, A}, pattern::BioSequence{A}; 
                            max_edits::Int64=2) where {S, A <: Alphabet}
    matches = approximate_search(fm, pattern; max_edits=max_edits)
    return [m.position for m in matches]
end

end


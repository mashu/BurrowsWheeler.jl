using BurrowsWheeler
using BioSequences
using Test

@testset "Approximate search" begin
    @testset "Basic approximate search" begin
        T = dna"ATCGATCG"
        fm = build_fmindex(T)
        
        # Exact match should be found with edit distance 0
        matches = approximate_search(fm, dna"ATCG"; max_edits=0)
        @test length(matches) == 2
        @test all(m.edit_distance == 0 for m in matches)
        @test matches[1].position == 1
        @test matches[2].position == 5
        
        # Pattern with one substitution (ATGG -> ATCG)
        matches = approximate_search(fm, dna"ATGG"; max_edits=1)
        @test length(matches) >= 1
        @test all(m.edit_distance <= 1 for m in matches)
        @test any(m.position == 1 && m.edit_distance == 1 for m in matches)
        
        # Pattern not in text, but with one edit allowed
        matches = approximate_search(fm, dna"AAAA"; max_edits=1)
        @test all(m.edit_distance <= 1 for m in matches)
    end
    
    @testset "Substitutions" begin
        T = dna"ATCGATCG"
        fm = build_fmindex(T)
        
        # Find ATCG with substitution in last position (ATCG -> ATCC)
        matches = approximate_search(fm, dna"ATCC"; max_edits=1)
        @test length(matches) >= 1
        @test any(m.edit_distance == 1 for m in matches)
        
        # Find ATCG with substitution in first position (ATCG -> CTCG)
        matches = approximate_search(fm, dna"CTCG"; max_edits=1)
        @test length(matches) >= 1
        @test any(m.edit_distance == 1 for m in matches)
    end
    
    @testset "Insertions and deletions" begin
        T = dna"ATCGATCG"
        fm = build_fmindex(T)
        
        # Pattern shorter than match (deletion in pattern)
        matches = approximate_search(fm, dna"ATC"; max_edits=1)
        @test length(matches) >= 1
        # Should find exact matches at positions 1, 5 and also with one deletion
        
        # Pattern longer (insertion in text)
        matches = approximate_search(fm, dna"ATCG"; max_edits=1)
        @test length(matches) >= 2  # Exact matches
        
        # Pattern with one base missing
        matches = approximate_search(fm, dna"TCG"; max_edits=1)
        @test length(matches) >= 1
    end
    
    @testset "Edit distance limits" begin
        T = dna"ATCGATCG"
        fm = build_fmindex(T)
        
        # With max_edits=0, should only find exact matches
        matches = approximate_search(fm, dna"ATGG"; max_edits=0)
        @test all(m.edit_distance == 0 for m in matches)
        exact_matches = approximate_search(fm, dna"ATCG"; max_edits=0)
        @test length(exact_matches) == 2
        @test all(m.edit_distance == 0 for m in exact_matches)
        
        # With max_edits=2, should find more matches
        matches_1 = approximate_search(fm, dna"ATGG"; max_edits=1)
        matches_2 = approximate_search(fm, dna"ATGG"; max_edits=2)
        @test length(matches_2) >= length(matches_1)
    end
    
    @testset "approximate_locate function" begin
        T = dna"ATCGATCG"
        fm = build_fmindex(T)
        
        # Should return positions
        positions = approximate_locate(fm, dna"ATCG"; max_edits=0)
        @test length(positions) == 2
        @test positions == [1, 5]
        
        # Should return positions for approximate matches
        positions = approximate_locate(fm, dna"ATGG"; max_edits=1)
        @test length(positions) >= 1
        @test all(p >= 1 && p <= length(T) for p in positions)
        
        # Should return positions for patterns that can match with edits
        positions = approximate_locate(fm, dna"TTTT"; max_edits=2)
        # Note: This might return some matches depending on implementation
        @test all(p >= 1 && p <= length(T) for p in positions)
    end
    
    @testset "RNA sequences" begin
        T = rna"AUAUAU"
        fm = build_fmindex(T)
        
        # Exact match
        matches = approximate_search(fm, rna"AUA"; max_edits=0)
        @test length(matches) == 2
        @test all(m.edit_distance == 0 for m in matches)
        
        # Approximate match with substitution
        matches = approximate_search(fm, rna"AAA"; max_edits=1)
        @test length(matches) >= 1
        @test all(m.edit_distance <= 1 for m in matches)
    end
    
    @testset "Edge cases" begin
        # Single character sequence
        T1 = dna"A"
        fm1 = build_fmindex(T1)
        matches = approximate_search(fm1, dna"A"; max_edits=0)
        @test length(matches) == 1
        @test matches[1].position == 1
        @test matches[1].edit_distance == 0
        
        # Pattern longer than sequence
        T2 = dna"AT"
        fm2 = build_fmindex(T2)
        matches = approximate_search(fm2, dna"ATCG"; max_edits=2)
        @test all(m.edit_distance <= 2 for m in matches)
        
        # Repeated characters
        T3 = dna"AAAA"
        fm3 = build_fmindex(T3)
        matches = approximate_search(fm3, dna"AAA"; max_edits=0)
        @test length(matches) == 2  # Positions 1 and 2
        @test all(m.edit_distance == 0 for m in matches)
    end
    
    @testset "Match correctness" begin
        T = dna"ATCGATCGATCG"
        fm = build_fmindex(T)
        
        # Verify that approximate matches are valid
        pattern = dna"ATGG"
        matches = approximate_search(fm, pattern; max_edits=1)
        
        for match in matches
            # Check that position is valid
            @test match.position >= 1
            @test match.position <= length(T)
            
            # Check that edit distance is within limit
            @test match.edit_distance <= 1
            
            # Check that alignment length is reasonable
            @test match.alignment_length >= 0
        end
        
        # Verify exact matches are found correctly
        exact_matches = approximate_search(fm, dna"ATCG"; max_edits=0)
        @test length(exact_matches) == 3
        @test [m.position for m in exact_matches] == [1, 5, 9]
        @test all(m.edit_distance == 0 for m in exact_matches)
    end
    
    @testset "Multiple edit operations" begin
        T = dna"ATCGATCG"
        fm = build_fmindex(T)
        
        # Pattern with 2 edits allowed
        matches = approximate_search(fm, dna"CCCC"; max_edits=2)
        # Should find matches if they exist within 2 edits
        @test all(m.edit_distance <= 2 for m in matches)
        
        # Pattern that needs 2 substitutions
        matches = approximate_search(fm, dna"TTTT"; max_edits=2)
        @test all(m.edit_distance <= 2 for m in matches)
    end
    
    @testset "Consistency with exact search" begin
        T = dna"ATCGATCGATCG"
        fm = build_fmindex(T)
        
        pattern = dna"ATCG"
        
        # Exact matches from approximate_search with max_edits=0 should match locate
        exact_positions = locate(fm, pattern)
        approx_matches = approximate_search(fm, pattern; max_edits=0)
        approx_positions = [m.position for m in approx_matches]
        
        @test sort(approx_positions) == sort(exact_positions)
        @test all(m.edit_distance == 0 for m in approx_matches)
    end
    
    @testset "Alignment length correctness" begin
        T = dna"ATCGATCG"
        fm = build_fmindex(T)
        
        # Exact match: alignment_length should equal pattern length
        matches = approximate_search(fm, dna"ATCG"; max_edits=0)
        for m in matches
            @test m.alignment_length == 4
            @test m.edit_distance == 0
        end
        
        # Substitution: alignment_length should equal pattern length
        matches = approximate_search(fm, dna"ATGG"; max_edits=1)
        for m in matches
            if m.edit_distance == 1
                @test m.alignment_length == 4  # 4 text chars consumed
            end
        end
        
        # Deletion: pattern shorter than text consumed
        # Pattern ATC (3 chars) matching ATCG (4 chars) with deletion
        matches = approximate_search(fm, dna"ATC"; max_edits=1)
        exact_matches = [m for m in matches if m.edit_distance == 0]
        if !isempty(exact_matches)
            @test exact_matches[1].alignment_length == 3  # 3 text chars consumed
        end
        
        # Insertion: pattern longer than text consumed
        # Pattern ATCCG (5 chars) matching ATCG (4 chars) with insertion
        matches = approximate_search(fm, dna"ATCCG"; max_edits=1)
        for m in matches
            if m.edit_distance == 1
                @test m.alignment_length == 4  # 4 text chars consumed, 5 pattern chars matched
            end
        end
    end
    
    @testset "Empty pattern" begin
        T = dna"ATCG"
        fm = build_fmindex(T)
        
        # Empty pattern should match all positions (with 0 edits)
        matches = approximate_search(fm, dna""; max_edits=0)
        # Empty pattern matches everywhere (but this might be implementation-specific)
        @test all(m.edit_distance == 0 for m in matches)
    end
    
    @testset "Combined edit operations" begin
        T = dna"ATCGATCG"
        fm = build_fmindex(T)
        
        # Pattern requiring substitution + insertion
        matches = approximate_search(fm, dna"ATCCG"; max_edits=2)
        @test all(m.edit_distance <= 2 for m in matches)
        
        # Pattern requiring substitution + deletion
        matches = approximate_search(fm, dna"ATCC"; max_edits=2)
        @test all(m.edit_distance <= 2 for m in matches)
        
        # Pattern requiring insertion + deletion
        matches = approximate_search(fm, dna"ATCCG"; max_edits=2)
        @test all(m.edit_distance <= 2 for m in matches)
        
        # Pattern requiring all three edit types
        matches = approximate_search(fm, dna"TTCCG"; max_edits=3)
        @test all(m.edit_distance <= 3 for m in matches)
    end
    
    @testset "Boundary cases" begin
        T = dna"ATCGATCG"
        fm = build_fmindex(T)
        
        # Pattern at start of text
        matches = approximate_search(fm, dna"ATC"; max_edits=1)
        @test any(m.position == 1 for m in matches)
        
        # Pattern at end of text
        matches = approximate_search(fm, dna"TCG"; max_edits=1)
        @test any(m.position == 6 for m in matches)  # Position 6 in "ATCGATCG"
        
        # Pattern spanning entire text
        matches = approximate_search(fm, dna"ATCGATCG"; max_edits=0)
        @test any(m.position == 1 && m.edit_distance == 0 for m in matches)
    end
    
    @testset "Duplicate position handling" begin
        T = dna"ATCGATCG"
        fm = build_fmindex(T)
        
        # A position should only appear once, with best (lowest) edit distance
        matches = approximate_search(fm, dna"ATCG"; max_edits=2)
        positions = [m.position for m in matches]
        
        # Check that each position appears only once
        @test length(positions) == length(unique(positions))
        
        # Check that best edit distance is kept for each position
        position_to_edit = Dict{Int64, Int64}()
        for m in matches
            if !haskey(position_to_edit, m.position) || m.edit_distance < position_to_edit[m.position]
                position_to_edit[m.position] = m.edit_distance
            end
        end
        
        for m in matches
            @test m.edit_distance == position_to_edit[m.position]
        end
    end
    
    @testset "Large edit distances" begin
        T = dna"ATCGATCGATCG"
        fm = build_fmindex(T)
        
        # Pattern requiring many edits (more than half pattern length)
        matches = approximate_search(fm, dna"TTTTTTTT"; max_edits=5)
        @test all(m.edit_distance <= 5 for m in matches)
        
        # max_edits equal to pattern length
        matches = approximate_search(fm, dna"ATCG"; max_edits=4)
        @test all(m.edit_distance <= 4 for m in matches)
        
        # max_edits larger than pattern length
        matches = approximate_search(fm, dna"AT"; max_edits=10)
        @test all(m.edit_distance <= 10 for m in matches)
    end
    
    @testset "Alignment length edge cases" begin
        T = dna"ATCGATCG"
        fm = build_fmindex(T)
        
        # Single character pattern
        matches = approximate_search(fm, dna"A"; max_edits=0)
        for m in matches
            @test m.alignment_length == 1
        end
        
        # Pattern with all deletions (pattern shorter than any match)
        matches = approximate_search(fm, dna"A"; max_edits=3)
        # Should find matches with various alignment lengths
        @test all(m.alignment_length >= 1 for m in matches)
        # Alignment length should be reasonable (not exceed text length from position)
        for m in matches
            # Position should be valid
            if m.position <= length(T)
                max_possible = length(T) - m.position + 1
                @test m.alignment_length <= max_possible
            end
            @test m.alignment_length >= 0
        end
        
        # Pattern with all insertions (pattern longer than match)
        matches = approximate_search(fm, dna"ATCGATCGATCG"; max_edits=4)
        # Alignment length should be reasonable (non-negative)
        for m in matches
            @test m.alignment_length >= 0
            @test m.edit_distance <= 4
        end
    end
    
    @testset "All DNA characters substitution" begin
        T = dna"ATCGATCG"
        fm = build_fmindex(T)
        
        # Test substitution with each DNA character
        for char in [DNA_A, DNA_T, DNA_C, DNA_G]
            pattern = LongSequence{DNAAlphabet{4}}([char, DNA_T, DNA_C, DNA_G])
            matches = approximate_search(fm, pattern; max_edits=1)
            # Should find matches with substitution
            @test all(m.edit_distance <= 1 for m in matches)
        end
    end
    
    @testset "Pattern matching entire text" begin
        T = dna"ATCG"
        fm = build_fmindex(T)
        
        # Pattern exactly matches text
        matches = approximate_search(fm, dna"ATCG"; max_edits=0)
        @test length(matches) == 1
        @test matches[1].position == 1
        @test matches[1].edit_distance == 0
        @test matches[1].alignment_length == 4
        
        # Pattern longer than text with edits
        matches = approximate_search(fm, dna"ATCGATCG"; max_edits=4)
        @test length(matches) >= 1
        @test any(m.position == 1 for m in matches)
    end
end


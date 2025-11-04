using BurrowsWheeler
using BioSequences
using Test

@testset "FM-index exact search" begin
    @testset "Basic search functionality" begin
        T = dna"ATAATA"
        fm = build_fmindex(T)
        
        # Test search function
        pattern1 = dna"ATA"
        sp1, ep1 = search(fm, pattern1)
        @test sp1 <= ep1
        @test count_occurrences(fm, pattern1) == 2
        positions1 = locate(fm, pattern1)
        @test length(positions1) == 2
        @test positions1 == [1, 4]
        
        pattern2 = dna"TA"
        sp2, ep2 = search(fm, pattern2)
        @test sp2 <= ep2
        @test count_occurrences(fm, pattern2) == 2
        positions2 = locate(fm, pattern2)
        @test length(positions2) == 2
        @test positions2 == [2, 5]
        
        pattern3 = dna"AA"
        sp3, ep3 = search(fm, pattern3)
        @test sp3 <= ep3
        @test count_occurrences(fm, pattern3) == 1
        positions3 = locate(fm, pattern3)
        @test length(positions3) == 1
        @test positions3 == [3]
        
        # Pattern not found
        pattern4 = dna"ATCG"
        sp4, ep4 = search(fm, pattern4)
        @test sp4 < 1 || sp4 > ep4
        @test count_occurrences(fm, pattern4) == 0
    end
    
    @testset "DNA sequences" begin
        T1 = dna"ATAATA"
        fm1 = build_fmindex(T1)
        
        @test count_occurrences(fm1, dna"ATA") == 2
        @test count_occurrences(fm1, dna"TA") == 2
        @test count_occurrences(fm1, dna"AA") == 1
        @test count_occurrences(fm1, dna"ATCG") == 0
        
        positions1 = locate(fm1, dna"ATA")
        @test positions1 == [1, 4]
        
        T2 = dna"ATCGATCG"
        fm2 = build_fmindex(T2)
        
        @test count_occurrences(fm2, dna"ATCG") == 2
        @test count_occurrences(fm2, dna"TCGA") == 1
        @test count_occurrences(fm2, dna"CGAT") == 1
        @test count_occurrences(fm2, dna"GATC") == 1
    end
    
    @testset "RNA sequences" begin
        T = rna"AUAUAU"
        fm = build_fmindex(T)
        
        @test count_occurrences(fm, rna"AUA") == 2
        @test count_occurrences(fm, rna"UAU") == 2
        @test count_occurrences(fm, rna"AU") == 3
        @test count_occurrences(fm, rna"ACGU") == 0
        
        positions = locate(fm, rna"AUA")
        @test length(positions) == 2
    end
    
    @testset "Edge cases" begin
        # Single character sequence
        T1 = dna"A"
        fm1 = build_fmindex(T1)
        @test count_occurrences(fm1, dna"A") == 1
        @test count_occurrences(fm1, dna"T") == 0
        
        # Repeated characters
        T2 = dna"AAAA"
        fm2 = build_fmindex(T2)
        @test count_occurrences(fm2, dna"A") == 4
        @test count_occurrences(fm2, dna"AA") == 3
        @test count_occurrences(fm2, dna"AAA") == 2
        @test count_occurrences(fm2, dna"AAAA") == 1
        
        # All suffixes
        T3 = dna"ATCG"
        fm3 = build_fmindex(T3)
        @test count_occurrences(fm3, dna"ATCG") == 1
        @test count_occurrences(fm3, dna"TCG") == 1
        @test count_occurrences(fm3, dna"CG") == 1
        @test count_occurrences(fm3, dna"G") == 1
    end
    
    @testset "Pattern matching correctness" begin
        T = dna"ATCGATCGATCG"
        fm = build_fmindex(T)
        
        pattern = dna"ATCG"
        sp, ep = search(fm, pattern)
        @test sp <= ep
        @test ep - sp + 1 == 3
        
        positions = locate(fm, pattern)
        @test length(positions) == 3
        @test positions == [1, 5, 9]
        
        # Verify that positions actually match the pattern
        for pos in positions
            @test T[pos:pos+3] == pattern
        end
    end
    
    @testset "Longer sequences" begin
        T = dna"ATCGATCGATCGATCGATCG"
        fm = build_fmindex(T)
        
        @test count_occurrences(fm, dna"ATCG") == 5
        @test count_occurrences(fm, dna"TCGA") == 4
        @test count_occurrences(fm, dna"CGAT") == 4
        @test count_occurrences(fm, dna"GATC") == 4
        
        positions = locate(fm, dna"ATCG")
        @test length(positions) == 5
        @test positions == [1, 5, 9, 13, 17]
    end
    
    @testset "Repeated patterns" begin
        T = dna"AAAAA"
        fm = build_fmindex(T)
        
        @test count_occurrences(fm, dna"A") == 5
        @test count_occurrences(fm, dna"AA") == 4
        @test count_occurrences(fm, dna"AAA") == 3
        @test count_occurrences(fm, dna"AAAA") == 2
        @test count_occurrences(fm, dna"AAAAA") == 1
    end
    
    @testset "Search range correctness" begin
        T = dna"ATAATA"
        fm = build_fmindex(T)
        
        pattern = dna"ATA"
        sp, ep = search(fm, pattern)
        @test sp >= 1
        @test ep >= sp
        @test ep <= length(fm.bwt)
        
        # Verify consistency between search, count_occurrences, and locate
        count = count_occurrences(fm, pattern)
        @test count == ep - sp + 1
        @test count == length(locate(fm, pattern))
    end
end

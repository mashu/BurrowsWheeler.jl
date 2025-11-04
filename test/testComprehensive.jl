using BurrowsWheeler
using BioSequences
using Test

@testset "Comprehensive FM-index tests" begin
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
        T1 = dna"A"
        fm1 = build_fmindex(T1)
        @test count_occurrences(fm1, dna"A") == 1
        @test count_occurrences(fm1, dna"T") == 0
        
        T2 = dna"AAAA"
        fm2 = build_fmindex(T2)
        @test count_occurrences(fm2, dna"A") == 4
        @test count_occurrences(fm2, dna"AA") == 3
        @test count_occurrences(fm2, dna"AAA") == 2
        @test count_occurrences(fm2, dna"AAAA") == 1
        
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
        
        count = count_occurrences(fm, pattern)
        @test count == ep - sp + 1
        @test count == length(locate(fm, pattern))
    end
end


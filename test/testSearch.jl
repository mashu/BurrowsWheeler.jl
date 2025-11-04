using BurrowsWheeler
using BioSequences
using Test

T = dna"ATAATA"
fm = build_fmindex(T)

@testset "FM-index search" begin
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
    
    pattern4 = dna"ATCG"
    sp4, ep4 = search(fm, pattern4)
    @test sp4 < 1 || sp4 > ep4
    @test count_occurrences(fm, pattern4) == 0
end


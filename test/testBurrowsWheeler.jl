using BioSequences
using Test
BWT = [DNA_A, DNA_T, DNA_T, DNA_A, DNA_Gap, DNA_A, DNA_A]
T = dna"ATAATA"
Tg = dna"ATC-A"
@test bwt_naïve(T) == BWT
@test_throws ArgumentError bwt_naïve(Tg)
@test sa(T) == [7,6,3,4,1,5,2]
@test_throws ArgumentError sa(Tg)
@test bwt(T) == BWT
@test_throws ArgumentError bwt(Tg)
@test bwt_naïve(T) == bwt(T)
@test rank(bwt(T)) == ([0, 0, 1, 1, 0, 2, 3], Dict(DNA_Gap=>1,DNA_A=>4,DNA_T=>2))
ranks, counts = rank(bwt(T))
col = Dict(DNA_Gap => (0,1), DNA_T => (5,7), DNA_A => (1,5))
@test bwmcol(counts) == col
@test revbwt(bwt(T)) == T

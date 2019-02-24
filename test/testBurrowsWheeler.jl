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

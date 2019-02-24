using BioSequences
using Test
BWT = [DNA_A, DNA_T, DNA_T, DNA_A, DNA_Gap, DNA_A, DNA_A]
@test bwt_naïve(dna"ATAATA") == BWT
@test_throws ArgumentError bwt_naïve(dna"ATC-A")
@test sa(dna"ATAATA") == [7,6,3,4,1,5,2]
@test_throws ArgumentError sa(dna"ATC-A")
@test bwt(dna"ATAATA") == BWT
@test_throws ArgumentError bwt(dna"ATC-A")
@test bwt_naïve(dna"ATAATA") == bwt(dna"ATAATA")

using BioSequences
using Test
BWT = [DNA_A, DNA_T, DNA_T, DNA_A, DNA_Gap, DNA_A, DNA_A]
@test bwt(dna"ATAATA") == BWT
@test_throws ErrorException bwt(dna"ATC-A")

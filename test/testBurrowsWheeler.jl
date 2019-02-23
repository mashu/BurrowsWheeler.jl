using BioSequences
using Test
T = dna"ATAATA"
BWT = [DNA_A, DNA_T, DNA_T, DNA_A, DNA_Gap, DNA_A, DNA_A]
@test bwt(T) == BWT

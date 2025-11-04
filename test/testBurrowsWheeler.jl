using BioSequences
using Test
using BurrowsWheeler

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

# Tests for _alphabet_type
@testset "_alphabet_type" begin
    @test BurrowsWheeler._alphabet_type(RNA) == RNAAlphabet{4}
    @test BurrowsWheeler._alphabet_type(DNA) == DNAAlphabet{4}
end

# Tests for revbwt with alphabet_type parameter
@testset "revbwt with alphabet_type" begin
    T_dna = dna"ATAATA"
    bw_dna = bwt(T_dna)
    @test revbwt(bw_dna, DNAAlphabet{4}) == T_dna
    
    T_rna = rna"AUAUAU"
    bw_rna = bwt(T_rna)
    @test revbwt(bw_rna, RNAAlphabet{4}) == T_rna
    
    # Test with different alphabet type explicitly
    @test revbwt(bw_dna, DNAAlphabet{4}) == T_dna
end

# Tests for _encode_sequence generic version
@testset "_encode_sequence generic" begin
    # Test with DNAAlphabet{4} (should use specific method, but verify it works)
    seq_dna = dna"ATCG"
    encoded_dna = BurrowsWheeler._encode_sequence(seq_dna)
    @test encoded_dna == [0x01, 0x04, 0x02, 0x03]  # A, T, C, G
    
    # Test with RNAAlphabet{4} (should use specific method, but verify it works)
    seq_rna = rna"AUCG"
    encoded_rna = BurrowsWheeler._encode_sequence(seq_rna)
    @test encoded_rna == [0x01, 0x04, 0x02, 0x03]  # A, U, C, G
    
    # Test generic version would need a different alphabet type
    # For now, we verify the DNA/RNA specific versions work correctly
end

# Tests for _encode_dna - all return branches
@testset "_encode_dna all branches" begin
    @test BurrowsWheeler._encode_dna(DNA_Gap) == 0x00
    @test BurrowsWheeler._encode_dna(DNA_A) == 0x01
    @test BurrowsWheeler._encode_dna(DNA_C) == 0x02
    @test BurrowsWheeler._encode_dna(DNA_G) == 0x03
    @test BurrowsWheeler._encode_dna(DNA_T) == 0x04
    
    # Test else branch (0xFF) - use an invalid DNA character
    # DNA_N or other ambiguous characters should trigger else branch
    # Note: DNA_N might not exist, so we'll test with a character that doesn't match any condition
    # Actually, DNA_N does exist in BioSequences, let's check if it triggers else
    # For now, we can test that any non-standard DNA character returns 0xFF
    # Since BioSequences might have DNA_N, let's test with that
    @test BurrowsWheeler._encode_dna(DNA_N) == 0xFF
end

# Tests for _encode_rna - all return branches
@testset "_encode_rna all branches" begin
    @test BurrowsWheeler._encode_rna(RNA_Gap) == 0x00
    @test BurrowsWheeler._encode_rna(RNA_A) == 0x01
    @test BurrowsWheeler._encode_rna(RNA_C) == 0x02
    @test BurrowsWheeler._encode_rna(RNA_G) == 0x03
    @test BurrowsWheeler._encode_rna(RNA_U) == 0x04
    
    # Test else branch (0xFF) - use an invalid RNA character
    # RNA_N or other ambiguous characters should trigger else branch
    @test BurrowsWheeler._encode_rna(RNA_N) == 0xFF
end

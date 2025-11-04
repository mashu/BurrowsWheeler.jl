# BurrowsWheeler.jl

High-performance Burrows-Wheeler Transform and FM-index implementation in Julia for genomic sequence analysis with exact and approximate pattern matching.

## Installation

```julia
using Pkg
Pkg.add("BurrowsWheeler")
```

## Quick Start

```julia
using BurrowsWheeler
using BioSequences

# Create a DNA sequence
sequence = dna"ATAATA"

# Build FM-index for fast pattern searching
fm = build_fmindex(sequence)

# Exact pattern matching
sp, ep = search(fm, dna"ATA")  # Returns (start_position, end_position)
count = count_occurrences(fm, dna"ATA")  # Returns 2
positions = locate(fm, dna"ATA")  # Returns [1, 4]

# Approximate pattern matching (with substitutions, insertions, deletions)
matches = approximate_search(fm, dna"ATGG"; max_edits=1)  # Find matches with up to 1 edit
best_positions = approximate_locate(fm, dna"ATGG"; max_edits=1)  # Returns [1, 5]
```

## Examples

### Burrows-Wheeler Transform

```julia
using BurrowsWheeler
using BioSequences

sequence = dna"ATAATA"

# Compute suffix array using O(n) SA-IS algorithm
sa_result = sa(sequence)

# Compute BWT from suffix array
bwt_result = bwt(sequence)

# Reverse BWT to recover original sequence
original = revbwt(bwt_result)  # Returns dna"ATAATA"
```

### Working with Different Sequence Types

```julia
using BurrowsWheeler
using BioSequences

# DNA sequences
dna_seq = dna"ATCGATCG"
fm_dna = build_fmindex(dna_seq)
count_occurrences(fm_dna, dna"ATCG")  # Returns 2

# RNA sequences
rna_seq = rna"AUAUAU"
fm_rna = build_fmindex(rna_seq)
count_occurrences(fm_rna, rna"AUA")  # Returns 2
```

### Approximate Pattern Matching

Find approximate matches allowing substitutions, insertions, and deletions:

```julia
using BurrowsWheeler
using BioSequences

# Build FM-index
genome = dna"ATCGATCGATCG"
fm = build_fmindex(genome)

# Search with up to 1 edit (substitution, insertion, or deletion)
pattern = dna"ATGG"
matches = approximate_search(fm, pattern; max_edits=1)

# Each match contains:
# - position: 1-based position in the original sequence
# - edit_distance: number of edits (substitutions, insertions, deletions)
# - alignment_length: length of the alignment

for match in matches
    println("Position: ", match.position, ", Edit distance: ", match.edit_distance)
end

# Or just get positions
positions = approximate_locate(fm, pattern; max_edits=1)
```

## API Reference

### Core Functions

```@autodocs
Modules = [BurrowsWheeler]
Order   = [:function, :type]
```

### Search Functions

```@autodocs
Modules = [BurrowsWheeler.Search]
Order   = [:type, :function]
```

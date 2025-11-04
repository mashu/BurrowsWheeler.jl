using Pkg
Pkg.develop(path=joinpath(@__DIR__, ".."))
using BurrowsWheeler
using BioSequences
using BenchmarkTools

const SUITE = BenchmarkGroup()

println("Running BurrowsWheeler.jl benchmarks...")

# Small sequences
SUITE["small"] = BenchmarkGroup()
T_small = dna"ATAATA"
fm_small = build_fmindex(T_small)
SUITE["small"]["bwt"] = @benchmarkable bwt($T_small)
SUITE["small"]["sa"] = @benchmarkable sa($T_small)
SUITE["small"]["build_fmindex"] = @benchmarkable build_fmindex($T_small)
SUITE["small"]["search"] = @benchmarkable search($fm_small, dna"ATA")
SUITE["small"]["count_occurrences"] = @benchmarkable count_occurrences($fm_small, dna"ATA")
SUITE["small"]["locate"] = @benchmarkable locate($fm_small, dna"ATA")

# Medium sequences
SUITE["medium"] = BenchmarkGroup()
T_medium = dna"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
fm_medium = build_fmindex(T_medium)
SUITE["medium"]["bwt"] = @benchmarkable bwt($T_medium)
SUITE["medium"]["sa"] = @benchmarkable sa($T_medium)
SUITE["medium"]["build_fmindex"] = @benchmarkable build_fmindex($T_medium)
SUITE["medium"]["search"] = @benchmarkable search($fm_medium, dna"ATCG")
SUITE["medium"]["count_occurrences"] = @benchmarkable count_occurrences($fm_medium, dna"ATCG")
SUITE["medium"]["locate"] = @benchmarkable locate($fm_medium, dna"ATCG")

# Large sequences
SUITE["large"] = BenchmarkGroup()
T_large = repeat(dna"ATCG", 100)
fm_large = build_fmindex(T_large)
SUITE["large"]["bwt"] = @benchmarkable bwt($T_large)
SUITE["large"]["sa"] = @benchmarkable sa($T_large)
SUITE["large"]["build_fmindex"] = @benchmarkable build_fmindex($T_large)
SUITE["large"]["search"] = @benchmarkable search($fm_large, dna"ATCG")
SUITE["large"]["count_occurrences"] = @benchmarkable count_occurrences($fm_large, dna"ATCG")
SUITE["large"]["locate"] = @benchmarkable locate($fm_large, dna"ATCG")

# Run benchmarks
results = run(SUITE, verbose=true, seconds=10)

println("\nBenchmark Results:")
println("==================")
println(results)

# Check if any benchmark is too slow (threshold: 10 seconds for large sequences)
for (size, group) in results
    for (op, trial) in group
        if median(trial).time > 10_000_000_000  # 10 seconds in nanoseconds
            println("⚠️  WARNING: $size/$op is slow: $(median(trial).time / 1e9) seconds")
        end
    end
end


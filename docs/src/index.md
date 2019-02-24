# BurrowsWheeler

Set of functions operating on BioSequence type from BioJulia.

```@docs
sa(word::BioSequence{DNAAlphabet{4}})::Array{Int64}
```

```@docs
bwt(word::BioSequence{DNAAlphabet{4}})::Array{DNA}
```

```@docs
bwt_naïve(word::BioSequence{DNAAlphabet{4}})::Array{DNA}
```

```@docs
rank(bw::Array{DNA})::Dict{DNA,Int64}
```

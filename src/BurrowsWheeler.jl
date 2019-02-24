module BurrowsWheeler
using BioSequences
export bwt, sa
    function sa(word::BioSequence{DNAAlphabet{4}})
        w = copy(word)
        if DNA_Gap in w
            throw(ArgumentError("input must not contain gaps"))
        end
        push!(w, DNA_Gap)
        S = sort!([(i, w[i:end]) for i in 1:(length(w))], by = x -> x[2])
        return(first.(S))
    end
    function bwt(word::BioSequence{DNAAlphabet{4}})
        if DNA_Gap in word
            throw(ArgumentError("input must not contain gaps"))
        end
        w = convert(Vector{DNA}, word)
        push!(w, DNA_Gap)
        l = Vector{DNA}()
        i = length(w)
        while i > 0
            c = pop!(w)
            pushfirst!(w, c)
            i = i - 1
            append!(l, w)
        end
        M = reshape(l, length(w), length(w))
        return(sortslices(M, dims=2)[end, :])
    end
end

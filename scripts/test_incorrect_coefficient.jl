using DrWatson
@quickactivate :MEK9470

##
function correct_coefficient(a_P, a_F)
    dP = 1 / a_P
    dF = 1 / a_F
    return 0.5(dP + dF)
end

function incorrect_coefficient(a_P, a_F)
    a_f = 0.5(a_P + a_F)
    return 1 / a_f
end

## d is diag of top left block of system matrix
correct_coefs = Float64[]
incorrect_coefs = Float64[]
for i in eachindex(d), j in eachindex(d)
    if i != j
        push!(correct_coefs, correct_coefficient(d[i], d[j]))
        push!(incorrect_coefs, incorrect_coefficient(d[i], d[j]))
    end
end

extrema(correct_coefs)
extrema(incorrect_coefs)

diffs = correct_coefs - incorrect_coefs
extrema(diffs)

reldiffs = diffs ./ correct_coefs
extrema(reldiffs)
function second_derivative(f::Vector{Float64}, h::Float64)
    n = length(f)
    d2f = zeros(n)
    for i in 2:n-1
        d2f[i] = (f[i-1] - 2 * f[i] + f[i+1]) / h^2
    end
    # Optional: Set boundary values to NaN or extrapolate
    d2f[1] = NaN
    d2f[end] = NaN
    return d2f
end
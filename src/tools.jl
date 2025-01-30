function create_1D_grid(A, B, N)
    dx = (B - A)/N
    x = collect(LinRange(A +  dx/2, B - dx/2, N))
    return x, dx
end


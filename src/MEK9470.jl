module MEK9470
    using LinearAlgebra
    using CairoMakie
    using Statistics
    export create_1D_grid
    export AbstractDiscretization, CentralDiff, Upwind, Hybrid, QUICK, TVD
    export VanLeer, VanAlbada, Minmod, Superbee, Sweby, UMIST, UD
    export apply_limiter
    include("tools.jl")
    include("types.jl")
end

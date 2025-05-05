module MEK9470
    using Reexport
    using LinearAlgebra
    using CairoMakie
    using Statistics


    export create_1D_grid
    include("tools.jl")
    include("Discretizations.jl/Discretizations.jl")
    @reexport using .Discretizations

    using NonlinearSolve
    using LinearAlgebra

    #LDCSolver
    using NonlinearSolve
    using DifferentialEquations
    using LinearAlgebra
    using LoopVectorization
    using LinearSolve
    using SparseArrays
    
    include("LDCSolver/types.jl")
    export AbstractMesh, CartesianMesh, AbstractLDCProblem, LDCProblem
    export LDCOperators, LinearOperator, ConvectionOperator
    include("LDCSolver/mesh_utils.jl")
    export linear_to_cartesian_indices, cartesian_to_linear_indices, cell_center
    include("LDCSolver/solver.jl")
    export laplacian, get_laplacian_operator, get_continuity_operator, continuity_x, continuity_y, pressure_dx, pressure_dy
    export get_pressure_dx_operator, get_pressure_dy_operator
    
end

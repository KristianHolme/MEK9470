module LDCSolver
    using NonlinearSolve
    using LinearAlgebra


    export AbstractMesh, CartesianMesh, AbstractLDCProblem, LDCProblem
    export linear_to_cartesian_indices

    include("types.jl")
    include("mesh_utils.jl")
    include("solver.jl")
    include("discretizations.jl")
    include()

end
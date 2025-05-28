module MEK9470

using Reexport
using LinearAlgebra
using SparseArrays
using CairoMakie
using Statistics
using NonlinearSolve
using DifferentialEquations
using LoopVectorization
using LinearSolve

# Basic tools and utilities
export create_1D_grid
include("tools.jl")

# Discretizations submodule
include("Discretizations/Discretizations.jl")
@reexport using .Discretizations

# Mesh types and utilities
include("LDCSolver/types.jl")
export AbstractMesh, CartesianMesh, AbstractLDCProblem, LDCProblem
export LDCOperators, LinearOperator, ConvectionOperator

include("LDCSolver/mesh_utils.jl")
export linear_to_cartesian_indices, cartesian_to_linear_indices, cell_center
export get_face_index, get_total_faces, face_neighbors

# Solver functionality
include("LDCSolver/solver.jl")
export laplacian, get_laplacian_operator, get_continuity_operator
export continuity_x, continuity_y, pressure_dx, pressure_dy
export get_pressure_dx_operator, get_pressure_dy_operator
export compute_face_values, compute_r_ratio, get_convection_operators

# Boundary condition system
export BoundaryCondition, DirichletBC, NeumannBC
export DomainBoundaryConditions, North, South, East, West
export lid_driven_cavity_u_bc, lid_driven_cavity_v_bc
export get_bc, get_boundary_side, get_boundary_value

# Face and flow types
export FaceOrientation, Vertical, Horizontal
export FlowDirection, Positive, Negative
export upwind_cell, downwind_cell, flow_direction, face_velocity

end

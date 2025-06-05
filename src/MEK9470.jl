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
using JLD2
using Dates
using SparseConnectivityTracer, ADTypes
using NonlinearSolve.NonlinearSolveBase: RelNormSafeBestTerminationMode

# Basic tools and utilities
export create_1D_grid
include("tools.jl")


include("LDCSolver/mesh_utils.jl")
export linear_to_cartesian_indices, cartesian_to_linear_indices, cell_center
export get_face_index, get_total_faces, face_neighbors
# Boundary condition system
export BoundaryCondition, DirichletBC, NeumannBC
export DomainBoundaryConditions, North, South, East, West
export lid_driven_cavity_u_bc, lid_driven_cavity_v_bc
export get_bc, get_boundary_side, get_boundary_value

# Face and flow types
export FaceOrientation, Vertical, Horizontal
export FlowDirection, Positive, Negative
export upwind_cell, downwind_cell, flow_direction_sign, face_velocity


# Discretizations submodule
include("Discretizations/Discretizations.jl")
@reexport using .Discretizations

# Solver functionality
include("LDCSolver/operators.jl")
export set_convection_operators!, update_convection_operators!
export laplacian, get_laplacian_operator, get_continuity_operator_x, get_continuity_operator_y
export continuity_x, continuity_y, pressure_dx, pressure_dy
export get_pressure_dx_operator, get_pressure_dy_operator
export compute_face_values, compute_r_ratio, get_convection_operators
export split_uvp

# Mesh types and utilities
include("LDCSolver/types.jl")
export AbstractMesh, CartesianMesh, AbstractLDCProblem, LDCProblem
export LDCOperators, LinearOperator, ConvectionOperator



include("LDCSolver/solver.jl")
export construct_linear_system, residuals

include("LDCPlotting/plotting.jl")
export simple_uvp_plot, lid_driven_cavity_plot, plot_velocity_profiles, plot_streamfunction, plot_flow_analysis

# Flow analysis utilities
include("LDCUtils.jl/flow_functions.jl")
export cumsimp, compute_streamfunction, compute_velocity_potential, compute_flow_functions
export compute_circulation, compute_vorticity, verify_flow_relations, compute_stream

include("LDCUtils.jl/solving_utils.jl")
export compute_next_re, solve_up

end

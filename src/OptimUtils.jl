module OptimUtils

# BuildUp tracker to log optimization summaries
using Base
using AbstractTrees
include("buildup.jl")
export BuildUp, NoBuildUp, OptionalBuildUp
export @addnode, @addbranch, addnode, sumall!, headnode, nodetype, innertree
include("inertia.jl")
export InertiaBuildUp, InertialElement

# Variable and output name, bounds and indices utilities
using StaticArrays
using Printf
include("var.jl")
export Var, ini, lower, upper, indexbyname, indexbygroup, len, mergevar, ini_scaled, get_scaled
export unscale_unpack, unscale_unpack!, views, unpack, unpack_s, getvar!, scale, subset, setvariables!

# Problem definition and convenient autodiff
using ForwardDiff
include("problem.jl")
end

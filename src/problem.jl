abstract type AbstractProblem end

struct ProblemCache{
    Prb<:AbstractProblem, 
    var<:NamedTuple, out<:NamedTuple,
    IDX<:NamedTuple, IDG<:NamedTuple,
    C,CD,opt,
    M<:DiffResults.MutableDiffResult,
    J<:ForwardDiff.JacobianConfig, T
    }
    # Problem
    prb::Prb

    # Variables
    variables::var
    x::Vector{T}
    x_lb::Vector{T}
    x_ub::Vector{T}
    idx::IDX

    # Objective and constraints
    output::out
    g::Vector{T}
    idg::IDG

    # Other inputs to analysis function
    usercache::C
    options::opt
    
    # Differentiation tools
    usercache_dual::CD
    fdresults::M
    cfg::J

    function ProblemCache(prb::Prb, variables::var, output::out, options::opt, chunksize::Int=len(variables)) where {Prb<:AbstractProblem,var,out,opt}
        # Var basics
        idx = indexbyname(variables)
        x = ini_scaled(variables)
        x_lb = lower(variables)
        x_ub = upper(variables)
        idg = indexbyname(output)
        g = ini(output)
        usercache = makecache(eltype(x), prb, options)
    
        # Dual 
        fdresults = DiffResults.JacobianResult(g, x)
        cfg = ForwardDiff.JacobianConfig(OptimUtils.analysis, g, x, ForwardDiff.Chunk{chunksize}(), nothing)
        usercache_dual = makecache(eltype(cfg), prb, options)
        
        new{Prb,var,out,typeof(idx),typeof(idg),typeof(usercache),typeof(usercache_dual),opt,typeof(fdresults),typeof(cfg),eltype(x)}(prb, variables, x, x_lb, x_ub, idx, output, g, idg, 
                            usercache, options, usercache_dual, fdresults, cfg)
    end
end

Base.show(io::IO, p::ProblemCache{Prb}) where Prb = print(io,
"""
ProblemCache for $Prb
Variables: $(keys(p.variables))
Output: $(keys(p.output))
""")
Base.show(io::IO, p::Type{<:ProblemCache{Prb}}) where Prb = print(io,"ProblemCache{$Prb}")


function analysis(p::ProblemCache, g::AbstractVector, x::AbstractVector, memory=nothing)
    throw("analysis function must be defined")
end

function makecache(T::DataType, p::AbstractProblem, options)
    nothing
end

function optfun(p::ProblemCache, g::AbstractVector, df::AbstractVector, dg::AbstractMatrix, x::AbstractVector, objective::Symbol)
    # possible improvement: reinterpret cache instead of disabling tag checking
    ForwardDiff.jacobian!(p.fdresults, (g, x) -> analysis(p, g, x), g, x, p.cfg, Val{false}())
    g .= p.fdresults.value
    dg .= p.fdresults.derivs[1]

    f = p.fdresults.value[p.idg[objective]]
    for i=eachindex(df)
        df[i] = p.fdresults.derivs[1][p.idg[objective], i]
    end
    return f
end

function optfun(p::ProblemCache, g::AbstractVector, dg::AbstractMatrix, x::AbstractVector)
    # possible improvement: reinterpret cache instead of disabling tag checking
    ForwardDiff.jacobian!(p.fdresults, (g, x) -> analysis(p, g, x), g, x, p.cfg, Val{false}())
    g .= p.fdresults.value
    dg .= p.fdresults.derivs[1]
    nothing
end
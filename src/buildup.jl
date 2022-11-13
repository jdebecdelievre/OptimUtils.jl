"""
Simple tree structure that allows to build drag and mass buildups easily.
Main constructor:
    tree = BuildUp(name, value)
Supports:
    - get child: child = tree[name]
    - add node from name, value: tree[name] = value
    - macros for efficient node or branch addition
        - @addnode tree node1 node2 node3...
        - @addbranch tree node subnode1 subnode2...
    - compute all node values from leaf information: sumall!. 
    - get headnode (cuts off children to avoid recomputing sums with all leaves): headnode(tree)
    - get single branch : branch(tree, new_head_node_name)
    - pretty printing print_tree
The node's value field can have any type. If using custom type, make sure to define sum() and + operations. 
See InertiaBuildUp for an example of custom value type.
"""

# Main node type:
mutable struct BuildUp{T}
    name::Symbol
    value::T
    children::Dict{Symbol,BuildUp{T}}
    @inline function BuildUp(name::Symbol, value::T) where T
        return new{T}(name, value, Dict{Symbol,BuildUp{T}}())
    end
end
Base.eltype(::BuildUp{T}) where T = T
Base.getindex(tree::BuildUp{T} where T, c::Symbol) = tree.children[c]
Base.setindex!(tree::BuildUp{T}, child::BuildUp{T}, name::Symbol) where T = setindex!(tree.children, child, name)
Base.setindex!(tree::BuildUp{T}, child::T, name::Symbol) where T = setindex!(tree.children, BuildUp(name,child), name)
Base.setindex!(tree::BuildUp{T}, child::Tv, name::Symbol) where {T,Tv} = setindex!(tree.children, BuildUp(name,T(child)), name)
branch(tree::BuildUp, c::Symbol) = tree.children[c]
headnode(tree::BuildUp) = BuildUp(tree.name, tree.value) # cuts off depth
function addnode(tree::BuildUp{T}, name::Symbol, value::Tv=tree.value; index::Union{Nothing,Int64}=nothing) where {T,Tv}
    name = isa(index, Nothing) ? name : Symbol("$(name)$index")
    (T!=Tv) && (value=T(value))
    tree.children[name]=BuildUp(name, value)
end

function Base.copy(tree::BuildUp)
    tnew = BuildUp(tree.name, tree.value)
    children = keys(tree.children) 
    for child in children
        tnew[child] = copy(tree[child])
    end
    return tnew
end

# Pretty printing. AbstractTrees only used for printing, dependency could be removed if necessary
AbstractTrees.printnode(io::IO, node::BuildUp) = print(io, "$(node.name): $(node.value)")
AbstractTrees.children(itree::BuildUp) = values(itree.children)
AbstractTrees.nodetype(::BuildUp) = BuildUp
Base.show(io::IO, tree::BuildUp) = print_tree(io, tree, depth=3) # using depth of 3

# Reduction operations
function reduce!(op, tree::BuildUp{T}, root::Bool=false) where T
    root && (tree.value = zero(T))
    if ~isempty(tree.children)
        if length(tree.children)==1
            child = iterate(values(tree.children))[1]
            tree.value = reduce!(op, child)
        else
            tree.value = op((reduce!(op, node) for node=values(tree.children))...)
        end
    end
    return tree.value
end
sumall!(tree::BuildUp) = reduce!(+,tree,true)
maxall!(tree::BuildUp) = reduce!(max,tree,true)
minall!(tree::BuildUp) = reduce!(min,tree,true)

"""
Binary operations on build-ups (+, -, max, min).
Both trees must be strictly identical: structure, names, and orders in which keys were added (that's due to the underlying use of dictionaries).
"""
function oper(op, t1::BuildUp{T}, t2::BuildUp{T}) where T
    tnew = BuildUp(t1.name, op(t1.value, t2.value))
    children = keys(t1.children) 
    if children == keys(t2.children)
        for child in children
            tnew[child] = oper(op, t1[child], t2[child])
        end
    end
    return tnew
end
Base.:+(t1::BuildUp{T}, t2::BuildUp{T})  where T = oper(+, t1, t2)
Base.:-(t1::BuildUp{T}, t2::BuildUp{T})  where T = oper(-, t1, t2)
Base.:*(t1::BuildUp{T}, t2::BuildUp{T})  where T = oper(*, t1, t2)
Base.:/(t1::BuildUp{T}, t2::BuildUp{T})  where T = oper(/, t1, t2)
Base.max(t1::BuildUp{T}, t2::BuildUp{T}) where T = oper(max, t1, t2)
Base.min(t1::BuildUp{T}, t2::BuildUp{T}) where T = oper(min, t1, t2)

function  Base.:^(tree::BuildUp{T}, n::Integer) where T 
    if n<1
        raise("Error: nonpositive power undefined for BuildUp trees")
    elseif n==1
        return tree
    else
        return tree * (tree)^(n-1)
    end
end

"""
Apply similar operation to all elements
"""
function scalar_oper!(op!, t::BuildUp{T} where T)
    t.value = op!(t.value)
    for tc=values(t.children)
        scalar_oper!(op!,tc)
    end
    return t
end
scalar_oper(op, tree::BuildUp{T} where T) = scalar_oper!(op, copy(tree))
# (op::Function)(tree::BuildUp{T} where T) = scalar_oper!(op, copy(tree)) # not working

"""
Binary operation between BuildUp{T} and object type T, or a Number (e.g multiplication by a scalar).
Note that the underlying type T in BuildUp{T} must have defined a the operation between T and T, T and a scalar.
"""
Base.:*(α::Union{T,Number}, tree::BuildUp{T}) where T = scalar_oper(t->α*t, tree)
Base.:*(tree::BuildUp{T}, α::Union{T,Number}) where T = scalar_oper(t->t*α, tree)
Base.:/(tree::BuildUp{T}  where T, α::Number) = tree*(1/α)
Base.:/(α::Number, tree::BuildUp{T}  where T) = scalar_oper(t->α/t, tree)

"""
Options for skipping buildup.
Example function signature:
    function f(inputs, buildup::OptionalBuildUp{T}=NoBuildUp)
        ...
    end
T is the type of the "value" field in your BuildUp, which could be a custom type (see InertialElement for an example)
"""
const NoBuildUp = BuildUp(:NoBuildUp, nothing)
const OptionalBuildUp{T} = Union{BuildUp{Nothing},BuildUp{T}}
Base.show(io::IO, tree::BuildUp{Nothing}) = Base.show(io,nothing)

# Duplicate all methods to make sure logging is ignored if necessary
Base.getindex(tn::BuildUp{<:Nothing}, ::Symbol) = nothing
Base.setindex!(tn::BuildUp{<:Nothing}, ::Any, ::Symbol) = tn
addnode(tn::BuildUp{<:Nothing}, ::Symbol, ::Any; index::Union{Nothing,Int64}=nothing) = tn
addnode(tn::BuildUp{T}, ::Symbol, ::T; index::Union{Nothing,Int64}=nothing) where T<:Nothing = tn
branch(::BuildUp{<:Nothing}, ::Symbol) = nothing
headnode(tn::BuildUp{<:Nothing}) = tn
reduce!(op,tn::BuildUp{<:Nothing}, root::Bool=false) = nothing
oper!(::BuildUp{T} where T,::BuildUp{<:Nothing}) = nothing
oper!(::BuildUp{<:Nothing}, ::BuildUp{T} where T) = nothing

## Macros for convenience
"""
```julia_skip
vnode = @addnode tree v1 v2 v3
```
Add one or more children to the tree whose names are v1,v2,..., and value are the values of v1,v2,....
Return newly created node, or nothing if tree was of type Nothing.
"""
macro addnode(tree, V...)
    T = Expr(:block, Tuple(:(OptimUtils.addnode($tree, $(Expr(:quote, v)), $v)) for v=V)...)
    esc(T)
end

"""
```julia_skip
V = @addbranch tree V(v1,v2,v3)
```
Add a child to the tree node Node whose name is V, and children are v1, v2, v3. 
Also creates nodes for v1,v2,v3 with respective names and values. 
Returns root node of subtree, or nothing if tree was of type Nothing.
"""
macro addbranch(tree, V)
    esc(quote
        vr = OptimUtils.addnode($tree, $(V.args)[1], ($tree).value)
        @addnode vr $(V.args[2:end]...)
    end)
end

"""
```julia_skip
innertree(tree, :mass)
```
When a custom type T was created for the "value" field of Buildup{T}, where T have several fields (such as mass), this function
returns the tree for only one of the dimensions (for instance mass).
See InertiaBuildUp for an example of custom T.
"""
function innertree(tree::BuildUp{T} where T, property::Symbol)
    node = BuildUp(tree.name, getproperty(tree.value, property)) 
    for c = values(tree.children)
        node[c.name] = innertree(c, property)
    end
    node
end
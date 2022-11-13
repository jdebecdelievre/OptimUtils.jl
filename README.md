# OptimUtils.jl
Helper tools for optimization.
Two main tools: BuildUp and Var
This package is constantly changing.
# BuildUp

This simple utility package exports the `BuildUp` type, which helps keep track of the source of an integrated quantity.

## Use cases

- track time spent on the various legs of a race when simulating a car
- track time spent on various calculations
- track sources of mass, inertia and center of gravity offset when designing an object using multiple elements
- track sources of drag when computing the total drag of an airplane

## Synthax

Main constructor:
- tree = BuildUp(name, value)
- OR a shortcut 
- tree = BuildUp(name, value)
  
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

# Var

Var allows to unpack variable vector into namedtuples of variables.
It also allows to clearly keep track of bounds (active and inactive).
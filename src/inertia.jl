using StaticArrays
"""
Inertial element. The matrix of inertia stored fields are with respect to the center of mass.
"""
struct InertialElement{TF}
    mass::TF
    xcg::TF
    ycg::TF
    zcg::TF
    Ixx::TF
    Iyy::TF
    Izz::TF
    Ixy::TF
    Iyz::TF
    Izx::TF
end
InertialElement(mass=0., xcg=zero(mass), ycg=zero(mass), zcg=zero(mass), 
                         Ixx=zero(mass), Iyy=zero(mass), Izz=zero(mass),
                         Ixy=zero(mass), Iyz=zero(mass), Izx=zero(mass)) = InertialElement(promote(mass, xcg, ycg, zcg, Ixx, Iyy, Izz, Ixy, Iyz, Izx)...)

Base.copy(IE::InertialElement) = InertialElement(IE.mass, IE.xcg, IE.ycg, IE.zcg, IE.Ixx, IE.Iyy, IE.Izz, IE.Ixy, IE.Iyz, IE.Izx)
Base.zero(::Union{Type{InertialElement{T}},InertialElement{T}}) where {T} = InertialElement(zero(T))
Base.show(io::IO, IE::InertialElement) = Base.show(io, tuple(IE.mass, IE.xcg, IE.ycg, IE.zcg, IE.Ixx, IE.Iyy, IE.Izz, IE.Ixy, IE.Iyz, IE.Izx))

"""
Moments of inertia with respect to any point
"""
function momentsofinertia(IR::InertialElement, x::Number=IE.xcg, y::Number=IE.ycg, z::Number=IE.zcg)
    Ixx = IE.Ixx + IE.mass * ((IE.ycg-y)^2 + (IE.zcg-z)^2)
    Iyy = IE.Iyy + IE.mass * ((IE.zcg-z)^2 + (IE.xcg-x)^2)
    Izz = IE.Izz + IE.mass * ((IE.ycg-y)^2 + (IE.xcg-x)^2)
    Ixy = IE.Ixy - IE.mass * (IE.xcg-y)*(IE.ycg-y)
    Iyz = IE.Iyz - IE.mass * (IE.ycg-y)*(IE.zcg-z)
    Izx = IE.Izx - IE.mass * (IE.zcg-z)*(IE.xcg-x)
    return Ixx, Iyy, Izz, Ixy, Iyz, Izx
end
function inertiamatrix(IR::InertialElement, x::Number=IE.xcg, y::Number=IE.ycg, z::Number=IE.zcg) 
    Ixx, Iyy, Izz, Ixy, Iyz, Izx = momentsofinertia(IR, x, y, z)
    return SA[Ixx Ixy Izx
              Ixy Iyy Iyz
              Izx Iyz Izz] 
end

function Base.:(+)(a::InertialElement{T1} where T1, b::InertialElement{T2} where T2)
    mass = a.mass + b.mass
    
    xcg = (a.xcg*a.mass + b.xcg*b.mass)/mass
    ycg = (a.ycg*a.mass + b.ycg*b.mass)/mass
    zcg = (a.zcg*a.mass + b.zcg*b.mass)/mass

    Ixx = (a.Ixx + a.mass * ((a.ycg-ycg)^2 + (a.zcg-zcg)^2)
         + b.Ixx + b.mass * ((b.ycg-ycg)^2 + (b.zcg-zcg)^2))
        
    Iyy = (a.Iyy + a.mass * ((a.zcg-zcg)^2 + (a.xcg-xcg)^2)
         + b.Iyy + b.mass * ((b.zcg-zcg)^2 + (b.xcg-xcg)^2))

    Izz = (a.Izz + a.mass * ((a.ycg-ycg)^2 + (a.xcg-xcg)^2)
         + b.Izz + b.mass * ((b.ycg-ycg)^2 + (b.xcg-xcg)^2))

    Ixy = (a.Ixy - a.mass * (a.xcg-ycg)*(a.ycg-ycg)
         + b.Ixy - b.mass * (b.xcg-ycg)*(b.ycg-ycg))
        
    Iyz = (a.Iyz - a.mass * (a.ycg-ycg)*(a.zcg-zcg)
         + b.Iyz - b.mass * (b.ycg-ycg)*(b.zcg-zcg))

    Izx = (a.Izx - a.mass * (a.zcg-zcg)*(a.xcg-xcg)
         + b.Izx - b.mass * (b.zcg-zcg)*(b.xcg-xcg))

    return InertialElement(mass, xcg, ycg, zcg, Ixx, Iyy, Izz, Ixy, Iyz, Izx)
end
InertiaBuildUp(name::Symbol) = BuildUp(name, InertialElement())

"""
Translate item by Δx,Δy or Δz.
(Ixx, Iyy, Izz, Ixy, Iyz, Izx) unchanged because it remains around the c.g.
"""
function translate(IE::InertialElement; Δx=0., Δy=0., Δz=0.)
    return InertialElement(IE.mass, IE.xcg+Δx, IE.ycg+Δy, IE.zcg+Δz, IE.Ixx, IE.Iyy, IE.Izz, IE.Ixy, IE.Iyz, IE.Izx)
end

"""
Rotate item around center of gravity by matrix R.
(Ixx, Iyy, Izz, Ixy, Iyz, Izx) rotated around the CG.
WARNING: This function needs testing
"""
function rotate(IE::InertialElement, R::AbstractMatrix)
    xcg, ycg, zcg = R * SA[IE.xcg, IE.ycg, IE.zcg]
    I = SA[IE.Ixx IE.Ixy IE.Izx
           IE.Ixy IE.Iyy IE.Iyz
           IE.Izx IE.Iyz IE.Izz]
    I= R*I*R'
    Ixx, Ixy, Izx = I[SA[1,2,3]]
    Iyy, Iyz = I[SA[5,6]]
    Izz = I[9]
    return InertialElement(IE.mass, xcg, ycg, zcg, Ixx, Iyy, Izz, Ixy, Iyz, Izx)
end

"""
Axis angle rotation. Axis in degrees
"""
function rotate(IE::InertialElement, a::AbstractVector, θ::Number)
    na = sqrt(sum(aa^2 for aa=a))
    K = SA[
        0.    -a[3]  a[2]
         a[3]    0. -a[1]
        -a[2]  a[1]    0.
    ]
    K = K / na
    Id = SA[
        1. 0. 0.
        0. 1. 0.
        0. 0. 1.
    ]
    R = Id + sind(θ) * K + (1-cosd(θ)) * (K*K)
    return rotate(IE, R)
end
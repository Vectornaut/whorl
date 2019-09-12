module PoincaréDisk

using Compose

export
  möbius_map,
  pts_to_pts,
  stable,
  planeproj,
  geodesic,
  ideal_edges,
  ideal_path,
  horotriangle,
  geodesic_orbit

# === möbius transformations

##[type cleanup] switch to fixed arrays!

# apply a möbius transformation, given as an operator on C^2, to a point on the
# complex plane
möbius_map(m::Matrix{T}, z) where T <: Number =
  (m[1,1]*z + m[1,2]) / (m[2,1]*z + m[2,2])

# find the derivative of a möbius transformation, given as an operator on C^2,
# at a point on the complex plane
function möbius_deriv(m::Matrix{T}, z) where T <: Number
  u = m[2,1]*z + m[2, 2]
  (m[1,1]*m[2,2] - m[1,2]*m[2,1]) / (u * u)
end

# get the möbius transformation
#  0 --> a
#  1 --> b
#  ∞ --> c
std_to_pts(a, b, c) = [[c*(b - a), b - a] [a*(c - b), c - b]]

# get the möbius transformation
#  a0 --> a1
#  b0 --> b1
#  c0 --> c1
pts_to_pts(a0, b0, c0, a1, b1, c1) =
  std_to_pts(a1, b1, c1) * inv(std_to_pts(a0, b0, c0))

# find the stable line of an element of GL(2,C)
function stable(m::Matrix{T}) where T <: Number
  # when you pass a matrix of type Hermitian to eigfact!, it calls LAPACK's
  # sygvd function, which puts the eigenvalues in ascending order. that means
  # the first eigenvector is the one that shrinks the most.
  ##eigvecs(Hermitian(m' * m))[:,1]
  eig(Hermitian(m' * m))[2][:,1] # eigvecs for Hermitian matrices seems to be broken now...
end

# === points, geodesics, and horocycles

# project a line in affine space to the complex plane
planeproj(v::Vector{T}) where T <: Number = v[1] / v[2]

# a version of reim for use in Compose paths
reim_measure(z::Number) = (real(z)*cx, imag(z)*cy)

# write the curveto command for a geodesic between two points on the boundary of
# the Poincaré disk
function geodesic_curveto(tail::Number, head::Number)
  # this is Clinton Curry's method for approximating geodesics by cubic curves
  # http://clintoncurry.nfshost.com/math/poincare-geodesics.html
  # the _max_ in k is a kludge to avoid square roots of negative numbers
  k = 4/3 * (1/(1 + sqrt(max(1 - abs2(head + tail)/4, 0))) - 1/4)
  [:C, reim_measure(k*tail)..., reim_measure(k*head)..., reim_measure(head)...]
end

# draw the geodesic between two points on the boundary of the Poincaré disk
geodesic(tail::Number, head::Number) =
  compose(
    context(units=UnitBox(-1, -1, 2, 2)),
    path([
      :M; reim_measure(tail)...;
      geodesic_curveto(tail, head)...
    ])
  )

# draw an ideal polygon as a sequence of geodesics, good for stroking
function ideal_edges(verts::Number...)
  n = length(verts)
  cyc = i -> mod(i, n) + 1
  compose([geodesic(verts[cyc(i)], verts[cyc(i+1)]) for i in 0:(n-1)]...)
end

# draw an ideal polygon as a path, good for filling
function ideal_path(verts::Number...)
  n = length(verts)
  cyc = i -> mod(i, n) + 1
  compose(
    context(units=UnitBox(-1, -1, 2, 2)),
    path([
      :M; reim_measure(verts[1])...;
      [geodesic_curveto(verts[cyc(i)], verts[cyc(i+1)]) for i in 0:(n-1)]...
    ])
  )
end

# draw the circular arc from tail to head with the specified direction at the
# the tail. the length of the direction vector doesn't matter. unstable at large
# radii and near right angles, but our usage, in the function horoarc, avoids
# both conditions.
function arc(tail::Number, head::Number, dir::Number)
  # this is Aleksas Riškus's method for approximating circular arcs by cubic
  # curves, found in the paper "Approximation of a cubic Bezier curve by
  # circular arcs and vice versa".
  
  # the angle from the ray tail --> head to the tangent ray, given as a point on
  # the unit circle
  sweep = dir / (head - tail)
  sweep /= abs(sweep)
  
  # the magic number
  k = 4/3*(sqrt(2) - 1) * abs(imag(sweep)/real(sweep))
  
  # the radius of the arc
  r = abs((head - tail) / 2imag(sweep))
  
  # unit tangent vectors at the head and the tail
  tail_tan = dir / abs(dir)
  head_tan = tail_tan / (sweep * sweep)
  
  compose(
    context(units=UnitBox(-1, -1, 2, 2)),
    curve(reim(tail), reim(tail + k*r*tail_tan), reim(head - k*r*head_tan), reim(head))
  )
end

# draw an arc of a horocycle at osc, starting at the geodesic osc -- a and
# ending at the geodesic osc -- b, at distance _height_ from the edge of the
# contact triangle
function horoarc(osc::Number, a::Number, b::Number, height::Number, eps::Number = 1e-2)
  # find the tail of the desired arc on the standard triangle ∞, 0, 1
  std_tail = 1im * exp(height)
  
  # apply the mobius transformation
  #  ∞ --> osc
  #  0 --> a
  #  1 --> b
  # to get the desired arc on the triangle osc, a, b
  m = std_to_pts(a, b, osc)
  tail = möbius_map(m, std_tail)
  head = möbius_map(m, std_tail + 1)
  if abs2(head - tail) > eps*eps
    return arc(tail, head, möbius_deriv(m, std_tail))
  else
    return compose(context())
  end
end

# foliate the the osc corner of the triangle osc, a, b using _cnt_ evenly spaced
# horocycles
horoleaves(osc::Number, a::Number, b::Number, cnt::Integer, spacing::Number, eps::Number = 1e-2) =
  compose(context(), [horoarc(osc, a, b, n*spacing, eps) for n in 0:(cnt - 1)]...)

# draw a triangle foliated by horocycles
horotriangle(osc::Number, a::Number, b::Number, cnt::Integer, spacing::Number, eps::Number = 1e-2) =
  compose(
    context(),
    horoleaves(osc, a, b, cnt, spacing, eps),
    horoleaves(a, b, osc, cnt, spacing, eps),
    horoleaves(b, osc, a, cnt, spacing, eps),
  )

end

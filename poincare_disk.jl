module PoincaréDisk

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

using LinearAlgebra, Compose

# === möbius transformations

##[type cleanup] switch to fixed arrays!

# apply a möbius transformation, given as an operator on C^2, to a point on the
# complex plane
möbius_map(m::Union{AbstractMatrix, UniformScaling}, z) =
  (m[1,1]*z + m[1,2]) / (m[2,1]*z + m[2,2])

# find the derivative of a möbius transformation, given as an operator on C^2,
# at a point on the complex plane
function möbius_deriv(m::Union{AbstractMatrix, UniformScaling}, z)
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
function stable(m::Union{AbstractMatrix, UniformScaling})
  # by default, `eigen` sorts the eigenvalues lexicographically by their real
  # and imaginary parts. that means the first eigenvector of a Hermitian matrix
  # is the one that shrinks the most.
  eigen(Hermitian(m' * m)).vectors[:,1]
end

# === points, geodesics, and horocycles

# project a line in affine space to the complex plane
planeproj(v::Vector{<:Number}) = v[1] / v[2]

# the control point scale factor in Clinton Curry's method for approximating
# geodesics by cubic curves
#   https://web.archive.org/web/20161020210738/http://clintoncurry.nfshost.com/math/poincare-geodesics.html
# the _max_ in k is a kludge to avoid square roots of negative numbers
curryfactor(tail::Number, head::Number) =
  4/3 * (1/(1 + sqrt(max(1 - abs2(head + tail)/4, 0))) - 1/4)

# list the control points and endpoint of a geodesic between two points on the
# boundary of the Poincaré disk. the `Compose.bezigon` method takes a list of
# such lists as its `sides` argument
function geodesic_side(tail::Number, head::Number)
  k = curryfactor(tail, head)
  [reim(k*tail), reim(k*head), reim(head)]
end

# draw the geodesic between two points on the boundary of the Poincaré disk
function geodesic(tail::Number, head::Number)
  k = curryfactor(tail, head)
  compose(
    context(units=UnitBox(-1, -1, 2, 2)),
    curve(reim(tail), reim(k*tail), reim(k*head), reim(head))
  )
end

# arguments can be passed in arrays in order to perform multiple drawing operations
function geodesic(tails::Vector{<:Number}, heads::Vector{<:Number})
  ks = [curryfactor(tail, head) for (tail, head) in zip(tails, heads)]
  compose(
    context(units=UnitBox(-1, -1, 2, 2)),
    curve(reim.(tails), reim.(ks.*tails), reim.(ks.*heads), reim.(heads))
  )
end

# draw an ideal polygon as a sequence of geodesics, good for stroking
ideal_edges(verts::Vector{<:Number}) = geodesic(verts, circshift(verts, -1))

# arguments can be passed in arrays in order to perform multiple drawing operations
function ideal_edges(polyverts::Vector{<:Vector{<:Number}})
  tails = vcat(polyverts...)
  heads = vcat([circshift(verts, -1) for verts in polyverts]...)
  geodesic(tails, heads)
end

# draw an ideal polygon as a path, good for filling
ideal_path(verts::Vector{<:Number}) =
  compose(
    context(units=UnitBox(-1, -1, 2, 2)),
    bezigon(
      reim(first(verts)),
      [geodesic_side(tail, head) for (tail, head) in zip(verts, circshift(verts, -1))]
    )
  )

# arguments can be passed in arrays in order to perform multiple drawing operations
ideal_path(polyverts::Vector{<:Vector{<:Number}}) =
  compose(
    context(units=UnitBox(-1, -1, 2, 2)),
    bezigon(
      reim.(first.(polyverts)),
      [
        [geodesic_side(tail, head) for (tail, head) in zip(verts, circshift(verts, -1))]
        for verts in polyverts
      ]
    )
  )

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

# arguments can be passed in arrays in order to perform multiple drawing operations
function arc(tails::Vector{<:Number}, heads::Vector{<:Number}, dirs::Vector{<:Number})
  # the angles from the ray tail --> head to the tangent ray
  sweeps = dirs ./ (heads .- tails)
  sweeps ./= abs.(sweeps)
  
  # the magic numbers
  ks = 4/3*(sqrt(2) - 1) * abs.(imag.(sweeps)./real.(sweeps))
  
  # the radii of the arcs
  rs = abs.((heads .- tails) ./ 2imag.(sweeps))
  
  # unit tangent vectors at the head and the tail
  tail_tans = dirs ./ abs.(dirs)
  head_tans = tail_tans ./ (sweeps .* sweeps)
  
  compose(
    context(units=UnitBox(-1, -1, 2, 2)),
    curve(reim.(tails), reim.(tails .+ ks.*rs.*tail_tans), reim.(heads .- ks.*rs.*head_tans), reim.(heads))
  )
end

# foliate the the osc corner of the triangle verts = [osc, a, b] using _cnt_
# evenly spaced horocycles, starting with the edge of the contact triangle
horoleaves(verts::Vector{<:Number}, cnt::Integer, spacing::Number, eps::Number = 1e-2) =
  horoleaves([verts], cnt, spacing, eps)

# arguments can be passed in arrays in order to perform multiple drawing operations
function horoleaves(polyverts::Vector{<:Vector{<:Number}}, cnt::Integer, spacing::Number, eps::Number = 1e-2)
  # find the tails of the desired arcs on the standard triangle ∞, 0, 1
  std_tails = [1im * exp(n*spacing) for n in 0:(cnt - 1)]
  
  tails = nothing
  heads = nothing
  dirs = nothing
  for (osc, a, b) in polyverts
    # apply the mobius transformation
    #  ∞ --> osc
    #  0 --> a
    #  1 --> b
    # to get the desired arcs on the triangle osc, a, b
    m = std_to_pts(a, b, osc)
    for n in 1:cnt
      tail = möbius_map(m, std_tails[n])
      head = möbius_map(m, std_tails[n] + 1)
      if abs2(head - tail) > eps*eps
        if tails == nothing
          tails = [tail]
          heads = [head]
          dirs = [möbius_deriv(m, std_tails[n])]
        else
          push!(tails, tail)
          push!(heads, head)
          push!(dirs, möbius_deriv(m, std_tails[n]))
        end
      else
        break
      end
    end
  end
  arc(tails, heads, dirs)
end

# draw a triangle foliated by horocycles
horotriangle(verts::Vector{<:Number}, cnt::Integer, spacing::Number, eps::Number = 1e-2) =
  horoleaves([verts, circshift(verts, -1), circshift(verts, -2)], cnt, spacing, eps)

# arguments can be passed in arrays in order to perform multiple drawing operations
horotriangle(polyverts::Vector{<:Vector{<:Number}}, cnt::Integer, spacing::Number, eps::Number = 1e-2) =
  horoleaves(
    vcat([[verts, circshift(verts, -1), circshift(verts, -2)] for verts in polyverts]...),
    cnt, spacing, eps
  )

end

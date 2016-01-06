using Compose

# === möbius transformations

# apply a möbius transformation, given as an operator on C^2, to a point on the
# complex plane
möbius_map{T <: Number}(m::Matrix{T}, z) =
  (m[1,1]*z + m[1,2]) / (m[2,1]*z + m[2,2])

# find the derivative of a möbius transformation, given as an operator on C^2,
# at a point on the complex plane
function möbius_deriv{T <: Number}(m::Matrix{T}, z)
  u = m[2,1]*z + m[2, 2]
  (m[1,1]*m[2,2] - m[1,2]*m[2,1]) / (u * u)
end

# === geodesics and horocycles

# find the repelling fixed point of a translation of the Poincaré disk
function repeller{T <: Number}(m::Matrix{T})
  # i'm assuming eigfact returns the eigenvalues of a hermitian matrix in
  # ascending order. that seems to be true, but i can't see why: even for
  # hermitian matrices, eigfact calls LAPACK's geevx rather than heevx, so
  # there shouldn't be any guarantee on the eigenvalue ordering.
  line = eigvecs(Hermitian(m' * m))[:, 1]
  
  # we're safe from small denominators here, because the ratio is on the unit
  # circle
  line[1] / line[2]
end

# draw the geodesic between two points on the boundary of the Poincaré disk
function geodesic(tail::Number, head::Number)
  # this is Clinton Curry's method for approximating geodesics by cubic curves
  # http://clintoncurry.nfshost.com/math/poincare-geodesics.html
  # the _max_ in k is a kludge to avoid square roots of negative numbers
  k = 4/3 * (1/(1 + sqrt(max(1 - abs2(head + tail)/4, 0))) - 1/4)
  compose(
    context(units=UnitBox(-1, -1, 2, 2)),
    curve(reim(tail), reim(k*tail), reim(k*head), reim(head))
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
# ending at the geodesic osc -- b. as the height parameter varies from 1 to 0,
# the arc shrinks from an edge of the orthic triangle down to the point osc.
function horoarc(osc::Number, a::Number, b::Number, height::Number)
  # find head, tail, and direction of the desired arc on the standard triangle
  # 1, cis(π/3), cis(2π/3)
  std_dir = cis(π/2 + height*π/6)
  std_head = 1 + sqrt(3)*(1im + conj(std_dir))
  std_tail = 1 + sqrt(3)*(-1im + std_dir)
  
  # build the möbius transformation
  #  1         --> osc     1         -->  0 --> osc
  #  cis(2π/3) --> a    =  cis(2π/3) -->  1 --> a
  #  cis(4π/3) --> b       cis(4π/3) --> -1 --> b
  m = [[2a*b - osc*(a + b), a + b - 2osc] [osc*(b - a), b - a]] * [[-1im, sqrt(3)] [1im, sqrt(3)]]
  
  # apply the mobius transformation to get the desired arc on the triangle
  # osc, a, b
  arc(
    möbius_map(m, std_tail),
    möbius_map(m, std_head),
    möbius_deriv(m, std_tail) * std_dir
  )
end

# foliate the the osc corner of the triangle osc, a, b using the specified
# number of horocycles, _density_
horoleaves(osc::Number, a::Number, b::Number, density::Integer) =
  compose(context(), [horoarc(osc, a, b, h/density) for h in 1:density]...)

# draw a triangle foliated by horocycles
function horotriangle(
  osc::Number, a::Number, b::Number,
  density::Integer, sidecolor="black", horocolor="black"
)
  sides = compose(
    context(),
    geodesic(osc, a),
    geodesic(a, b),
    stroke(sidecolor)
  )
  bulk = compose(
    context(),
    horoleaves(osc, a, b, density),
    horoleaves(a, b, osc, density),
    horoleaves(b, osc, a, density),
    stroke(horocolor)
  )
  compose(context(), sides, bulk)
end

# === orbit drawing

# draw the orbit of a geodesic under a finitely generated group of isometries
function geodesic_orbit(
  tail::Number,
  head::Number,
  sym,
  depth::Integer,
  last_sym = nothing
)
  orbit = []
  
  # draw the geodesic
  push!(orbit, geodesic(tail, head))
  
  # if we're going deeper into the orbit, apply the generators and their
  # inverses and recurse
  if depth > 0
    for s in 1:length(sym)
      # apply generator
      if last_sym != -s
        f_head = möbius_map(sym[s], head)
        f_tail = möbius_map(sym[s], tail)
        push!(orbit, geodesic_orbit(f_tail, f_head, sym, depth - 1, s))
      end
      
      # apply inverse generator
      if last_sym != s
        b_head = möbius_map(inv(sym[s]), head)
        b_tail = möbius_map(inv(sym[s]), tail)
        push!(orbit, geodesic_orbit(b_tail, b_head, sym, depth - 1, -s))
      end
    end
  end
  
  compose(context(), orbit...)
end

# draw the orbit of a horocycle-foliated ideal triangle under a finitely
# generated group of isometries
function triangle_orbit(
  a::Number,
  b::Number,
  c::Number,
  density::Integer,
  sym,
  depth::Integer,
  last_sym = nothing
)
  orbit = []
  
  # draw the horocycle foliation
  push!(orbit, horoleaves(c, a, b, density))
  
  # draw the sides
  ##push!(orbit, compose(context(), stroke("orange"), geodesic(c, a)))
  push!(orbit, compose(context(), stroke("green"), geodesic(a, b)))
  push!(orbit, compose(context(), stroke("violet"), geodesic(b, c)))
  
  # if we're going deeper into the orbit, apply the generators and their
  # inverses and recurse
  if depth > 0
    for s in 1:length(sym)
      # apply generator
      if last_sym != -s
        fa = möbius_map(sym[s], a)
        fb = möbius_map(sym[s], b)
        fc = möbius_map(sym[s], c)
        push!(orbit, triangle_orbit(fa, fb, fc, density, sym, depth - 1, s))
      end
      
      # apply inverse generator
      if last_sym != s
        ba = möbius_map(inv(sym[s]), a)
        bb = möbius_map(inv(sym[s]), b)
        bc = möbius_map(inv(sym[s]), c)
        push!(orbit, triangle_orbit(ba, bb, bc, density, sym, depth - 1, -s))
      end
    end
  end
  
  compose(context(), orbit...)
end

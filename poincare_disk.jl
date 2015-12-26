using Compose

# apply a möbius transformation, given as an operator on C^2, to a complex
# number
möbius_map{T <: Number}(m::Matrix{T}, z) =
  (m[1,1]*z + m[1,2]) / (m[2,1]*z + m[2,2])

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
  # this is clinton curry's method for approximating geodesics by cubic curves
  # http://clintoncurry.nfshost.com/math/poincare-geodesics.html
  k = 4/3 * (1/(1 + sqrt(1 - abs2(head + tail)/4)) - 1/4)
  compose(
    context(units=UnitBox(-1, -1, 2, 2)),
    fill(nothing),
    curve(reim(tail), reim(k*tail), reim(k*head), reim(head))
  )
end

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

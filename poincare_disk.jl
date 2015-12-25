using Compose

# find the repelling fixed point of a translation of the Poincaré disk
function repeller{T <: Number}(tras::Matrix{T})
  # i'm assuming eigfact returns the eigenvalues of a hermitian matrix in
  # ascending order. that seems to be true, but i can't see why: even for
  # hermitian matrices, eigfact calls LAPACK's geevx rather than heevx, so
  # there shouldn't be any guarantee on the eigenvalue ordering.
  line = eigvecs(Hermitian(tras' * tras))[:, 1]
  
  # we're safe from small denominators here, because the ratio is on the unit
  # circle
  line[1] / line[2]
end

# draw the geodesic between two points on the boundary of the Poincaré disk
function geodesic(tail::Number, head::Number)
  arcangle = sqrt(head/tail)
  radius = abs(imag(arcangle)/real(arcangle))
  
  if real(arcangle) > 1e-6
    # the geodesic is visibly curved
    arc = path([
            :M, real(tail)*cx, imag(tail)*cy,
            :A, radius*cx, radius*cy, 0, false, imag(arcangle) < 0, real(head)*cx, imag(head)*cy
          ])
  else
    # the geodesic looks straight
    arc = line([(real(tail)*cx, imag(tail)*cy), (real(head)*cx, imag(head)*cy)])
  end
  
  compose(context(units=UnitBox(-1, -1, 2, 2)), fill(nothing), arc)
end

module Square

using ValidatedNumerics, IntervalExchange

export SquareLocSys, GenSquareLocSys, cocycle, appx_abelianize

abstract SquareLocSys

type GenSquareLocSys <: SquareLocSys
  n_transit
  e_transit
  s_transit
  w_transit
  
  GenSquareLocSys(n_transit, e_transit) = new(
    n_transit, e_transit,
    inv(n_transit), inv(e_transit)
  )
end

function cocycle{R <: AbstractInterval}(angle::R, loc::SquareLocSys)
  if (angle < @interval(pi)/2)
    off_perp = angle - @interval(pi)/4
    f_transit = [loc.n_transit, loc.e_transit]
  elseif (angle < @interval(pi))
    off_perp = angle - 3@interval(pi)/4
    f_transit = [loc.w_transit, loc.n_transit]
  elseif (angle < 3@interval(pi)/2)
    off_perp = angle - 5@interval(pi)/4
    f_transit = [loc.s_transit, loc.w_transit]
  else
    off_perp = angle - 7@interval(pi)/4
    f_transit = [loc.e_transit, loc.s_transit]
  end
  Cocycle([1 + tan(off_perp), 2], f_transit, [2, 1])
end

function abelianize{R <: AbstractInterval}(angle::R, loc::SquareLocSys, depth::Integer)
  # build and evolve cocycle
  orig = cocycle(angle, loc)
  iter = power_twostep(orig, depth)
  
  # abelianize
  ab = IntervalExchange.abelianize(orig, iter)
  ab_transit = [block.f_transit for block in ab.blocks_by_in]
  if (angle < @interval(pi)/2)
    GenSquareLocSys(ab_transit[1], ab_transit[2])
  elseif (angle < @interval(pi))
    GenSquareLocSys(ab_transit[2], inv(ab_transit[1]))
  elseif (angle < 3@interval(pi)/2)
    GenSquareLocSys(inv(ab_transit[1]), inv(ab_transit[2]))
  else
    GenSquareLocSys(inv(ab_transit[2]), ab_transit[1])
  end
end

# this approximation is accurate for angles very close to pi/2
function appx_abelianize(tilt_sgn, loc::SquareLocSys)
  # go into a basis where the north transition map is diagonal
  n_eigvals, n_eigframe = eig(loc.n_transit)
  if (abs2(n_eigvals[1]) < abs2(n_eigvals[2]))
    λ = n_eigvals[1]
  else
    λ = n_eigvals[2]
    n_eigframe = [n_eigframe[:,2] n_eigframe[:,1]]
  end
  e_reframed = inv(n_eigframe) * loc.e_transit * n_eigframe
  
  # approximate the east transition map of the abelianized local system
  if tilt_sgn > 0
    e_transit_ab = diagm([e_reframed[1,1], 1/e_reframed[1,1]])
  else
    e_transit_ab = diagm([1/e_reframed[2,2], e_reframed[2,2]])
  end
  
  # return
  GenSquareLocSys(diagm([λ, 1/λ]), e_transit_ab)
end

end

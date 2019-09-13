# this causes some kind of conflict when it's included twice, so you gotta do it
# by hand
##include("interval_exchange.jl")

module Rectangle

using ValidatedNumerics, IntervalExchange

export RectangleLocSys, cocycle, appx_abelianize

type RectangleLocSys{R <: AbstractInterval}
  width::R
  height::R
  n_transit
  e_transit
  s_transit
  w_transit
  
  RectangleLocSys(
    width, height,
    n_transit, e_transit
  ) = new(
    width, height,
    n_transit, e_transit,
    inv(n_transit), inv(e_transit)
  )
end

## Julia, i don't have time for your incomprehensible type system any more
#function cocycle{R <: AbstractInterval}(angle::R, loc::RectangleLocSys)
function cocycle(angle, loc)
  if (angle < @interval(pi)/2)
    measure = [loc.width * sin(angle), loc.height * cos(angle)]
    f_transit = [loc.e_transit * loc.n_transit, loc.e_transit]
  else
    ##measure = [-loc.height * cos(angle), loc.width * sin(angle)]
    ##f_transit = [loc.w_transit, loc.n_transit]
  end
  Cocycle(cumsum(measure), f_transit, [2, 1])
end

function abelianize{R <: AbstractInterval}(angle::R, loc::RectangleLocSys, depth::Integer)
  # build cocycle
  orig = cocycle(angle, loc)
  
  # evolve cocycle
  iter = orig
  for _ in 1:depth
    iter = twostep(iter)
  end
  
  # abelianize
  ab = IntervalExchange.abelianize(orig, iter)
  ab_transit = [block.f_transit for block in ab.blocks_by_in]
  if (angle < @interval(pi)/2)
    RectangleLocSys{AbstractInterval}(loc.width, loc.height, ab_transit[1], ab_transit[2])
  else
    RectangleLocSys{AbstractInterval}(loc.width, loc.height, ab_transit[2], inv(ab_transit[1]))
  end
end

# this approximation is accurate for angles very close to pi/2
function appx_abelianize{R <: AbstractInterval}(tilt_sgn, loc::RectangleLocSys{R})
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
    e_transit_ab = Diagonal([e_reframed[1,1], 1/e_reframed[1,1]])
  else
    e_transit_ab = Diagonal([1/e_reframed[2,2], e_reframed[2,2]])
  end
  
  # return
  RectangleLocSys{R}(
    loc.width,
    loc.height,
    Diagonal([λ, 1/λ]),
    e_transit_ab
  )
end

end

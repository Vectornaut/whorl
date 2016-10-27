# this causes some kind of conflict when it's included twice, so you gotta do it
# by hand
##include("interval_exchange.jl")

module Rectangle

using ValidatedNumerics, IntervalExchange

export RectangleLocSys, cocycle, appx_abelianize

type RectangleLocSys{R <: AbstractInterval}
  width::R
  height::R
  w_transit
  n_transit
  e_transit
  s_transit
  
  RectangleLocSys(
    width, height,
    w_transit, n_transit
  ) = new(
    width, height,
    w_transit, n_transit,
    inv(w_transit), inv(n_transit)
  )
end

function cocycle{R <: AbstractInterval}(angle::R, loc::RectangleLocSys{R})
  if (angle < @interval(pi)/2)
    measure = [loc.width * sin(angle), loc.height * cos(angle)]
    f_transit = [loc.n_transit, loc.e_transit]
  else
    measure = [-loc.height * cos(angle), loc.width * sin(angle)]
    f_transit = [loc.w_transit, loc.n_transit]
  end
  Cocycle(cumsum(measure), f_transit, [2, 1])
end

# this approximation is accurate for angles very close to pi/2
function appx_abelianize{R <: AbstractInterval}(tilt_sgn, loc::RectangleLocSys{R})
  # go into a basis where the north transition map is diagonal
  n_eigvals, n_eigframe = eig(loc.n_transit)
  if (abs2(n_eigvals[1]) < abs2(n_eigvals[2]))
    λ = n_eigvals[1]
  else
    λ = n_eigvals[2]
    n_eigframe = [n_eigframe[2] n_eigframe[1]]
  end
  w_reframed = inv(n_eigframe) * loc.w_transit * n_eigframe
  
  # approximate the west transition map of the abelianized local system
  if tilt_sgn > 0
    w_transit_ab = diagm([1/w_reframed[2,2], w_reframed[2,2]])
  else
    w_transit_ab = diagm([w_reframed[1,1], 1/w_reframed[1,1]])
  end
  
  # return
  RectangleLocSys{R}(
    loc.width,
    loc.height,
    w_transit_ab,
    diagm([λ, 1/λ])
  )
end

end

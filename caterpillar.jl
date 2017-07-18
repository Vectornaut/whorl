module Caterpillar

using ValidatedNumerics, IntervalExchange

export CaterpillarLocSys, cocycle, almost_flat_caterpillar

tilt(x::Integer, y::Integer) = ((10*x + 7*y)//149, (-7*x + 10*y)//149)

type CaterpillarLocSys
  f_transit
  
  CaterpillarLocSys(m4, m5, m7, m8, m9, a, b) = new(
    Array[
      b*inv(m9)*a,
      b*inv(m8)*a,
      b*inv(m7)*a,
      m9,
      m8,
      m7,
      inv(a)*b,
      m5,
      m4,
      inv(a)*inv(b),
      inv(a)*inv(m5)*inv(b),
      inv(a)*inv(m4)*inv(b),
      a*inv(b)
    ]
  )
end

function proj{R <: AbstractInterval, T <: Integer}(angle::R, pt::Tuple{T, T})
  # the rotation that lays the long diagonal of the caterpillar horizontal
  tilt = [[10, -7] [7, 10]]
  
  # the projection direction, in the tilted point of view
  tilt_up = tilt * [cos(angle), sin(angle)]
  
  # the point we want to project, in the tilted point of view
  tilt_pt = tilt * collect(pt)
  
  (tilt_pt[1] - (tilt_up[1]/tilt_up[2])*tilt_pt[2])::R
end

function cocycle{R <: AbstractInterval}(angle::R, loc::CaterpillarLocSys)
  # branch points of a presentation of a genus-five surface with four
  # singularities
  branch_pts = [(0, 1), (1, 1), (2, 2), (2, 3), (3, 3), (4, 4), (5, 4), (6, 5), (7, 5), (7, 6), (8, 7), (9, 7), (10, 7)]
  
  breaks = [proj(angle, pt) for pt in branch_pts]
  f_shuffle = [7, 12, 11, 10, 9, 8, 13, 6, 5, 4, 3, 2, 1]
  Cocycle(breaks, loc.f_transit, f_shuffle)
end

function almost_flat_caterpillar(g)
  m4 = g[4]
  m5 = g[3]
  m7 = g[2]
  m8 = g[1]
  m9 = m8*inv(m7)*m5*inv(m4)
  a = -m5*inv(m4)
  b = -inv(m4)*m5
  CaterpillarLocSys(m4, m5, m7, m8, m9, a, b)
end

end

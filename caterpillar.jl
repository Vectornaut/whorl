using ValidatedNumerics

include("interval_exchange.jl")
include("regular.jl")

tilt(x::Integer, y::Integer) = ((10*x + 7*y)//149, (-7*x + 10*y)//149)

function proj{T <: Integer}(angle::Interval, pt::Tuple{T, T})
  # the rotation that lays the long diagonal of the caterpillar horizontal
  tilt = [[10, -7] [7, 10]]
  
  # the projection direction, in the tilted point of view
  tilt_up = tilt * [cos(angle), sin(angle)]
  
  # the point we want to project, in the tilted point of view
  tilt_pt = tilt * collect(pt)
  
  tilt_pt[1] - (tilt_up[1]/tilt_up[2])*tilt_pt[2]
end

function caterpillar_cocycle(angle::Interval, m4, m5, m7, m8, m9, a, b)
  # branch points of a presentation of a genus-five surface with four
  # singularities
  branch_pts = [(0, 1), (1, 1), (2, 2), (2, 3), (3, 3), (4, 4), (5, 4), (6, 5), (7, 5), (7, 6), (8, 7), (9, 7), (10, 7)]
  breaks = Interval[proj(angle, pt) for pt in branch_pts]
  
  f_transit = Array[b*inv(m9)*a, b*inv(m8)*a, b*inv(m7)*a, m9, m8, m7, inv(a)*b, m5, m4, inv(a)*inv(b), inv(a)*inv(m5)*inv(b), inv(a)*inv(m4)*inv(b), a*inv(b)]
  f_shuffle = [7, 12, 11, 10, 9, 8, 13, 6, 5, 4, 3, 2, 1]
  Cocycle(breaks, f_transit, f_shuffle)
end

function twisted_caterpillar(angle::Interval)
  g = generators(2)
  m4 = g[4]
  m5 = g[3]
  m7 = g[2]
  m8 = g[1]
  m9 = m8*inv(m7)*m5*inv(m4)
  a = -m5*inv(m4)
  b = -inv(m4)*m5
  caterpillar_cocycle(angle, m4, m5, m7, m8, m9, a, b)
end

# === testing

function caterpillar_test(angle_change::Interval)
  a = twisted_caterpillar(@interval(3π/4) + angle_change)
  
  for h in a.blocks
    println(h)
  end
  
  for i in 1:3
    a = twostep(a)
  end
  
  draw(
    PDF("lam_test.pdf", 10cm, 10cm),
    compose(context(), lamination(a, generators(2), 1), stroke("black"), linewidth(0.1))
  )
end

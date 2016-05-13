include("interval_exchange.jl")
include("regular.jl")

module Testing

using ValidatedNumerics, Compose, Colors, IntervalExchange, Regular

tilt(x::Integer, y::Integer) = ((10*x + 7*y)//149, (-7*x + 10*y)//149)

function proj{R <: AbstractInterval, T <: Integer}(angle::R, pt::Tuple{T, T})
  # the rotation that lays the long diagonal of the caterpillar horizontal
  tilt = [[10, -7] [7, 10]]
  
  # the projection direction, in the tilted point of view
  tilt_up = tilt * [cos(angle), sin(angle)]
  
  # the point we want to project, in the tilted point of view
  tilt_pt = tilt * collect(pt)
  
  (tilt_pt[1] - (tilt_up[1]/tilt_up[2])*tilt_pt[2])::R
end

function caterpillar_cocycle{R <: AbstractInterval}(angle::R, m4, m5, m7, m8, m9, a, b)
  # branch points of a presentation of a genus-five surface with four
  # singularities
  branch_pts = [(0, 1), (1, 1), (2, 2), (2, 3), (3, 3), (4, 4), (5, 4), (6, 5), (7, 5), (7, 6), (8, 7), (9, 7), (10, 7)]
  breaks = [proj(angle, pt) for pt in branch_pts]
  
  f_transit = Array[b*inv(m9)*a, b*inv(m8)*a, b*inv(m7)*a, m9, m8, m7, inv(a)*b, m5, m4, inv(a)*inv(b), inv(a)*inv(m5)*inv(b), inv(a)*inv(m4)*inv(b), a*inv(b)]
  f_shuffle = [7, 12, 11, 10, 9, 8, 13, 6, 5, 4, 3, 2, 1]
  Cocycle(breaks, f_transit, f_shuffle)
end

function twisted_caterpillar{R <: AbstractInterval}(angle::R)
  g = Regular.generators(2)
  m4 = g[4]
  m5 = g[3]
  m7 = g[2]
  m8 = g[1]
  m9 = m8*inv(m7)*m5*inv(m4)
  a = -m5*inv(m4)
  b = -inv(m4)*m5
  caterpillar_cocycle(angle, m4, m5, m7, m8, m9, a, b)
end

# === output

function caterpillar_pics{R <: AbstractInterval}(angle_offset::R = @interval(1/11); svg = false)
  # set up cocycle
  a = twisted_caterpillar(@interval(3π/4) + angle_offset)
  for h in a.blocks_by_in
    println(h)
  end
  println("by in ---------")
  for h in a.blocks_by_out
    println(h)
  end
  println("by out ---------")
  for i in 1:8
    a = @time(twostep(a))
    println("$(length(a.blocks_by_in)) blocks by in")
    println("$(length(a.blocks_by_out)) blocks by out")
    println("$i ~~~~~~~~~")
    #=
    Profile.print(maxdepth=8)
    println("---------")
    =#
  end
  
  #=
  # draw poincaré disk
  disk = compose(context(), circle(), stroke("black"), fill(nothing), linewidth(0.4mm))
  
  # draw lamination and foliation
  clay = RGB(161/255, 149/255, 126/255)
  #silt = RGB(204/255, 193/255, 174/255)
  amethyst = RGB(204/255, 125/255, 189/255)
  lam = compose(context(), lamination(a, Regular.generators(2), 3), stroke(clay), linewidth(0.1mm))
  fol = compose(context(), foliage(a), stroke(amethyst), linewidth(0.1))
  
  # print outputs
  if svg
    lam_file = SVG("laminated.svg", 7cm, 7cm)
    fol_file = SVG("foliated.svg", 7cm, 7cm)
  else
    lam_file = PDF("laminated.pdf", 7cm, 7cm)
    fol_file = PDF("foliated.pdf", 7cm, 7cm)
  end
  draw(lam_file, compose(context(), disk, lam))
  draw(fol_file, compose(context(), disk, lam, fol))
  =#
end

end

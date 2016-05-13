include("interval_exchange.jl")
include("regular.jl")
include("cayley_crawler.jl")

module Testing

using ValidatedNumerics, Compose, Colors, IntervalExchange, PoincaréDisk, Regular, Crawl

export caterpillar_pics, triangle_orbiter

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

function dot_orbiter(m)
  w = möbius_map(m, 0)
  loc = compose(
    context(units=UnitBox(-1, -1, 2, 2)),
    (context(real(w) - 0.01, imag(w) - 0.01, 0.02, 0.02), circle())
  )
  
end

triangle_orbiter(z1, z2, z3) =
  m -> begin
    w1 = möbius_map(m, z1)
    w2 = möbius_map(m, z2)
    w3 = möbius_map(m, z3)
    return compose(
      context(),
      (context(),
        geodesic(w1, w2),
        geodesic(w2, w3),
        geodesic(w3, w1),
        stroke("orangered")
      ),
      (horotriangle(w1, w2, w3, 60, 1/11), stroke("orange")),
      linewidth(0.1mm)
    )
  end

function caterpillar_pics{R <: AbstractInterval}(angle_offset::R = @interval(1/11); svg = false)
  # set up cocycle
  a = twisted_caterpillar(@interval(3π/4) + angle_offset)
  for h in a.blocks
    println(h)
  end
  
  # retrieve a break point
  tripod_break = a.blocks[6].out_left
  
  # evolve cocycle
  println("---------")
  for i in 1:3
    a = twostep(a)
  end
  for h in a.blocks
    println(h)
  end
  
  # find the vertices of the triangle corresponding to tripod_break
  f_break = nothing
  for h in a.blocks
    if strictprecedes(h.in_left, tripod_break) && strictprecedes(tripod_break, h.in_right)
      f_break = h
      break
    end
  end
  b_breaks = []
  for k in a.blocks
    if !missed_connection(f_break, k)
      push!(b_breaks, k)
    end
  end
  vertices = [repeller(f_break.f_transit); [repeller(k.b_transit) for k in b_breaks]]
  println(vertices)
  
  # enumerate symmetry group elements
  transit = generators(2)
  transit = [transit; [inv(t) for t in transit]]
  crawler = CayleyCrawler(4, 4, 2)
  find_home!(crawler, transit)
  
  # draw poincaré disk
  disk = compose(context(),
    (context(), circle(), fill("white"), stroke(nothing)),
    (context(), rectangle(), fill("gainsboro"), stroke(nothing))
  )
  
  # draw dots
  glass = RGBA(0.0, 0.8, 0.6, 0.2)
  dots = compose(context(), fan_out(dot_orbiter, crawler)..., fill(glass))
  
  # draw triangle lifts
  triangles = fan_out(triangle_orbiter(vertices[1], vertices[length(vertices)], vertices[2]), crawler)
  lam = compose(context(), triangles...)
  
  draw(PDF("tripod_test.pdf", 7cm, 7cm), compose(lam, disk))
  draw(PDF("crawler_test.pdf", 7cm, 7cm), compose(dots, disk))
end

end

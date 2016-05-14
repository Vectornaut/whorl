include("regular.jl")
include("cayley_crawler.jl")
include("caterpillar.jl")

module Testing

using ValidatedNumerics, Compose, Colors, IntervalExchange, Caterpillar, PoincaréDisk, Regular, Crawl

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
  a = twisted_caterpillar(@interval(3π/4) + angle_offset, Regular.generators(2))
  for h in a.blocks_by_in
    println(h)
  end
  println("by in ---------")
  for h in a.blocks_by_out
    println(h)
  end
  println("by out ---------")
  
  # retrieve a break point
  tripod_break = a.blocks_by_in[6].out_left
  
  # evolve cocycle
  for i in 1:3
    a = @time(twostep(a))
    println("$(length(a.blocks_by_in)) blocks by in")
    println("$(length(a.blocks_by_out)) blocks by out")
    println("$i ~~~~~~~~~")
  end
  
  # find the vertices of the triangle corresponding to tripod_break
  f_break = nothing
  for h in a.blocks_by_in
    if strictprecedes(h.in_left, tripod_break) && strictprecedes(tripod_break, h.in_right)
      f_break = h
      break
    end
  end
  b_breaks = []
  for k in a.blocks_by_in
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
  findhome!(crawler, transit)
  
  # draw poincaré disk
  disk = compose(context(),
    (context(), circle(), fill("white"), stroke(nothing)),
    (context(), rectangle(), fill("gainsboro"), stroke(nothing))
  )
  
  # draw dots
  glass = RGBA(0.0, 0.8, 0.6, 0.2)
  dots = compose(context(), mapcollect(dot_orbiter, crawler)..., fill(glass))
  
  # draw triangle lifts
  triangles = mapcollect(triangle_orbiter(vertices[1], vertices[length(vertices)], vertices[2]), crawler)
  lam = compose(context(), triangles...)
  
  draw(PDF("tripod_test.pdf", 7cm, 7cm), compose(lam, disk))
  draw(PDF("crawler_test.pdf", 7cm, 7cm), compose(dots, disk))
end

end

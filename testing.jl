include("regular.jl")
include("cayley_crawler.jl")
include("caterpillar.jl")

module Testing

using ValidatedNumerics, Compose, Colors, IntervalExchange, Caterpillar, PoincaréDisk, Regular, Crawl

type Triangle
  a::Complex
  b::Complex
  pivot::Complex
  gap::Real
  sing
  
  Triangle(a, b, pivot, sing) = new(a, b, pivot, 1-real(b/a), sing)
end

Base.isless(p::Triangle, q::Triangle) = p.gap < q.gap

function dot_orbiter(m)
  w = möbius_map(m, 0)
  compose(
    context(units=UnitBox(-1, -1, 2, 2)),
    (context(real(w) - 0.01, imag(w) - 0.01, 0.02, 0.02), circle())
  )
end

triangle_orbiter(p::Triangle) =
  m -> begin
    a = möbius_map(m, p.a)
    b = möbius_map(m, p.b)
    c = möbius_map(m, p.pivot)
    leafcolor = [
      "hotpink",
      "orangered",
      "gold",
      "purple"
    ][p.sing]
    return compose(
      context(),
      (context(),
        geodesic(a, b),
        geodesic(b, c),
        geodesic(c, a),
        stroke("gray")
      ),
      (horotriangle(a, b, c, 60, 1/11), stroke(leafcolor)),
      linewidth(0.1mm)
    )
  end

function test{R <: AbstractInterval}(angle_offset::R = @interval(1/11); svg = false)
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
  
  # hand-label singularities
  a.blocks_by_in[1].sing = 1
  a.blocks_by_in[2].sing = 2
  a.blocks_by_in[3].sing = 1
  a.blocks_by_in[4].sing = 2
  a.blocks_by_in[5].sing = 1
  a.blocks_by_in[6].sing = 2
  a.blocks_by_in[7].sing = 3
  a.blocks_by_in[8].sing = 4
  a.blocks_by_in[9].sing = 3
  a.blocks_by_in[10].sing = 4
  a.blocks_by_in[11].sing = 3
  a.blocks_by_in[12].sing = 4
  a.blocks_by_in[13].sing = 1
  
  # evolve cocycle
  for i in 1:3
    a = @time(twostep(a))
    println("$(length(a.blocks_by_in)) blocks by in")
    println("$(length(a.blocks_by_out)) blocks by out")
    println("$i ~~~~~~~~~")
  end
  
  # find the widest triangle for each singularity
  triangles = scancollect(a, Triangle,
    b_fn = (left, right, pivot) -> Triangle(
      repeller(right.b_transit),
      repeller(left.b_transit),
      repeller(pivot.f_transit),
      left.sing
    )
  )
  widest = Triangle[]
  for sing in 1:4
    push!(widest, maximum(filter(p -> p.sing == sing, triangles)))
  end
  
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
  
  #=
  # draw dots
  glass = RGBA(0.0, 0.8, 0.6, 0.2)
  dots = compose(context(), mapcollect(dot_orbiter, crawler)..., fill(glass))
  =#
  
  # draw triangle lifts
  tri = vcat([mapcollect(triangle_orbiter(p), crawler) for p in widest]...)
  tri_pic = compose(context(), tri...)
  
  #=
  # draw lamination and horocycle foliation
  lam = compose(context(), lamination(a, generators(2), 3)..., stroke("midnightblue"), linewidth(0.1mm))
  fol = compose(context(), foliage(a)...)
  =#
  
  draw(PDF("triangle_test.pdf", 7cm, 7cm), compose(tri_pic, disk))
  #=
  draw(PDF("crawler_test.pdf", 7cm, 7cm), compose(dots, disk))
  draw(PDF("laminated.pdf", 7cm, 7cm), compose(lam, disk))
  draw(PDF("foliated.pdf", 7cm, 7cm), compose(lam, fol, disk))
  =#
end

end

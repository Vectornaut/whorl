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
      RGB(255/255, 1/255, 73/255),
      RGB(255/255, 121/255, 1/255),
      RGB(255/255, 210/255, 0/255),
      #RGB(118/255, 200/255, 0/255)
      RGB(0/255, 200/255, 146/255)
      #RGB(98/255, 28/255, 158/255)
    ][p.sing]
    return compose(
      context(),
      (context(),
        ideal_path(a, b, c),
        fill(leafcolor)
      )
    )
  end

function test{R <: AbstractInterval}(angle_offset::R = @interval(1/11); frame = nothing, svg = false)
  # set up cocycle
  a = twisted_caterpillar(@interval(3π/4) + angle_offset, Regular.generators(2))
  #=
  for h in a.blocks_by_in
    println(h)
  end
  println("by in ---------")
  for h in a.blocks_by_out
    println(h)
  end
  println("by out ---------")
  =#
  
  # evolve cocycle
  for i in 1:3
    a = @time(twostep(a))
    println("$(length(a.blocks_by_in)) blocks")
    println("Step $i ~~~~~~~~~")
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
  
  # draw background
  bg = compose(context(),
    (context(), rectangle(), fill("white"), stroke(nothing))
  )
  
  # draw triangle lifts
  tri = vcat([mapcollect(triangle_orbiter(p), crawler) for p in widest]...)
  tri_pic = compose(context(), tri...)
  
  if frame == nothing
    draw(PDF("triangle_test.pdf", 7cm, 7cm), compose(tri_pic, bg))
  else
    draw(PNG(@sprintf("triangle_mov/frame%02i.png", frame), 500px, 500px), compose(tri_pic, bg))
  end
end

# the first two terms of a sawtooth wave, modified to zero out the jerk at the
# the inflection point and rescaled into the box [0,1] × [0,1]
function easing(t)
  θ = π*(t - 1/2)
  x = sin(θ) - sin(3θ)/26
  (1 + (x / (1+1/26))) / 2
end

function movie()
  start = @interval(-π/4 + 0.01)
  fin = @interval(π/4 - 0.01)
  n = 25
  for t in 0:n
    println(u)
    @time(test(u*start + (1-u)*fin, frame = t))
    println("Frame $t =========")
  end
end

end

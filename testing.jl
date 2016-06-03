include("regular.jl")
include("cayley_crawler.jl")
include("caterpillar.jl")

module Testing

using ValidatedNumerics, IntervalExchange, Caterpillar, PoincaréDisk, Regular, Crawl

type Triangle
  a::Complex
  b::Complex
  pivot::Complex
  gap::Real
  sing
  
  Triangle(a, b, pivot, sing) = new(a, b, pivot, 1-real(b/a), sing)
end

Base.isless(p::Triangle, q::Triangle) = p.gap < q.gap

function test()
  # set up cocycle
  a_orig = twisted_caterpillar(@interval(3π/4) + 1//11, Regular.generators(2))
  for b in a_orig.blocks_by_in
    println(b)
  end
  
  # evolve cocycle
  a = a_orig
  for i in 1:4
    a = @time(twostep(a))
    println("$(length(a.blocks_by_in)) blocks")
    println("$i ~~~~~~~~~")
  end
  
  println()
  a_ab = abelianize(a_orig, a)
  for i in 1:length(a_ab.blocks_by_in)
    println("block $i\n$(a_ab.blocks_by_in[i].f_transit)\n")
  end
end

end

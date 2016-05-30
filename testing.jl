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

type Jump{T <: Number, R <: AbstractInterval}
  op::Matrix{T}
  loc::R
  
  function Jump(left, right, pivot, loc)
    left_sc = left / det([left pivot])
    right_sc = right / det([right pivot])
    new([left_sc pivot] / [right_sc pivot], loc)
  end
end

function test{R <: AbstractInterval}(angle_offset::R = @interval(1/11))
  # set up cocycle
  a_orig = twisted_caterpillar(@interval(3π/4) + angle_offset, Regular.generators(2))
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
  
  # compute abelianization jumps
  ab_jumps = scancollect(a, Jump,
    f_fn = (left, right, pivot) -> Jump{Complex, Interval{Float64}}(
      stable(left.f_transit),
      stable(right.f_transit),
      stable(pivot.b_transit),
      left.in_right
    ),
    b_fn = (left, right, pivot) -> Jump{Complex, Interval{Float64}}(
      stable(left.b_transit),
      stable(right.b_transit),
      stable(pivot.f_transit),
      left.out_right
    )
  )
  
  println("\n$(length(ab_jumps)) jumps\n")
  
  big_bl = a_orig.blocks_by_in[4]
  y = big_bl.in_left + 0.01
  x = big_bl.out_left + 0.01
  println("x = $x\ny = $y\n")
  
  fine_bl_f = nothing
  fine_bl_b = nothing
  for b in a.blocks_by_in
    if strictprecedes(b.in_left, y) && strictprecedes(y, b.in_right)
      fine_bl_f = b
    end
    if strictprecedes(b.out_left, y) && strictprecedes(y, b.out_right)
      fine_bl_b = b
    end
    if fine_bl_f != nothing && fine_bl_b != nothing
      break
    end
  end
  println("fine blocks\n$fine_bl_f\n$fine_bl_b\ncontain y\n")
  println("forward-stable line: $(repeller(fine_bl_f.f_transit))")
  println("backward-stable line: $(repeller(fine_bl_b.b_transit))\n")
  
  dev_jumps = filter(
    j -> strictprecedes(y, j.loc) && strictprecedes(j.loc, x),
    ab_jumps
  )
  dev = prod([j.op for j in dev_jumps])
  hol = dev * big_bl.f_transit
  vals, vecs = eig(hol)
  println("eigenvalues:\n$(vals)")
  println("stable lines:\n$(vecs[1,1]/vecs[2,1])\n$(vecs[1,2]/vecs[2,2])")
  
  #=
  # enumerate symmetry group elements
  transit = generators(2)
  transit = [transit; [inv(t) for t in transit]]
  crawler = CayleyCrawler(4, 4, 2)
  findhome!(crawler, transit)
  =#
end

end

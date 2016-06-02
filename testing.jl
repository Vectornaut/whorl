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

type OldJump{T <: Number, R <: AbstractInterval}
  op::Matrix{T}
  loc::R
  
  function OldJump(left, right, pivot, loc)
    left_sc = left / det([left pivot])
    right_sc = right / det([right pivot])
    new([left_sc pivot] / [right_sc pivot], loc)
  end
end

function test(index::Integer)
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
  
  # compute abelianization jumps
  big_bl = a_orig.blocks_by_in[index]
  
  println("\n==== old code, block $index ====")
  ab_jumps = scancollect(a, OldJump,
    f_fn = (left, right, pivot) -> OldJump{Complex, Interval{Float64}}(
      stable(left.f_transit),
      stable(right.f_transit),
      stable(pivot.b_transit),
      left.in_right
    ),
    b_fn = (left, right, pivot) -> OldJump{Complex, Interval{Float64}}(
      stable(left.b_transit),
      stable(right.b_transit),
      stable(pivot.f_transit),
      left.out_right
    )
  )
  
  println("\n$(length(ab_jumps)) jumps\n")
  
  y = big_bl.in_left + 0.01
  x = big_bl.out_left + 0.01
  leftward = strictprecedes(y, x)
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
  
  jump_filter = leftward ?
    j -> strictprecedes(y, j.loc) && strictprecedes(j.loc, x) :
    j -> strictprecedes(x, j.loc) && strictprecedes(j.loc, y)
  dev_jumps = filter(jump_filter, ab_jumps)
  dev = prod([j.op for j in dev_jumps])
  if !leftward
    dev = inv(dev)
  end
  
  #=
  println("\n|||| new code, block $index ||||")
  ab_jumps = scancollect(a, Jump, f_fn = FJump, b_fn = BJump)
  dev_jumps = filter(
    j -> if isa(j, FJump)
      big_bl.orig_in <= j.left && j.pivot < big_bl.orig_out
    elseif isa(j, BJump)
      big_bl.orig_in <= j.pivot && j.left < big_bl.orig_out
    end,
    ab_jumps
  )
  dev = prod([j.op for j in dev_jumps])
  =#
  
  # output
  hol = dev * big_bl.f_transit
  vals, vecs = eig(hol)
  println("jumps: $(length(dev_jumps))")
  println("deviation:\n$dev\n")
  println("holonomy:\n$hol\n")
  println("eigenvalues:\n$vals\n")
  println("eigenlines:\n$(vecs[1,1]/vecs[2,1])\n$(vecs[1,2]/vecs[2,2])\n")
  
  println("\n++++ newer code, block $index ++++")
  a_ab = abelianize(a_orig, a, index)
  return
  #=
  for i in 1:length(a_ab.blocks_by_in)
    println("block $i\n$(a_ab.blocks_by_in[i].f_transit)\n")
  end
  =#
  
  #=
  # enumerate symmetry group elements
  transit = generators(2)
  transit = [transit; [inv(t) for t in transit]]
  crawler = CayleyCrawler(4, 4, 2)
  findhome!(crawler, transit)
  =#
end

end

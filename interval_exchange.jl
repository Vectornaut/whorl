using ValidatedNumerics, Compose

include("poincare_disk.jl")

# === exchangers

# an Exchanger is a piece of an interval exchange cocycle. it maps points in an
# "in" block to points in an "out" block by translation, discarding points that
# don't belong to the in block. it also applies a transition map to the
# local system sections above the block.
type Exchanger
  # the endpoints of the in and out blocks, specified using interval arithmetic.
  # the differences in_right - in_left and out_right - out_left will always be
  # equal, because the in block is mapped to the out block by translation.
  in_left::Interval
  in_right::Interval
  out_left::Interval
  out_right::Interval
  
  # the translation that sends the in block to the out block, always equal to
  # out_left - in_left
  f_shift::Interval
  
  # f_transit is the local system transition you pick up going forward through
  # the interval, and b_transit is the transition you pick up going backward.
  # f_transit and b_transit will always be inverses of one another.
  f_transit
  b_transit
  
  function Exchanger(in_left, in_right, f_shift, f_transit, b_transit)
    new(
      in_left,            # in_left
      in_right,           # in_right
      in_left + f_shift,  # out_left
      in_right + f_shift, # out_right
      f_shift,            # f_shift
      f_transit,          # f_transit
      b_transit,          # b_transit
    )
  end
end

# interval arithmetic collision exception
type EndpointCollisionException <: Exception
  h::Exchanger
  k::Exchanger
end

# checks whether the out block of k misses the in block of h
function missed_connection(h::Exchanger, k::Exchanger)
  # if we can't tell, throw an exception
  if !isdisjoint(h.in_right, k.out_left) || !isdisjoint(h.in_left, k.out_right)
    throw(EndpointCollisionException(h, k))
  end
  
  strictprecedes(h.in_right, k.out_left) || strictprecedes(k.out_right, h.in_left)
end

# compose the exchangers h and k, returning nothing if the out block of k
# misses the in block of h
function pipe(h::Exchanger, k::Exchanger)
  # if out block of k misses the in block of h, return nothing
  if missed_connection(h, k)
    return nothing
  end
  
  # find the left endpoint of the in block of the composition
  new_left = k.in_left
  if precedes(k.out_left, h.in_left) || !isdisjoint(k.out_left, h.in_left)
    # notice that this line runs if the left endpoint of the in block of h
    # collides with the left endpoint of the out block of k. that's because the
    # collision contributes to the uncertainty in the left endpoint of the in
    # block of the composition
    new_left += h.in_left - k.out_left
  end
  
  # find the right endpoint of the in block of the composed exchanger
  new_right = k.in_right
  if precedes(h.in_right, k.out_right) || !isdisjoint(h.in_right, k.out_right)
    # notice that this line runs if the right endpoint of the in block of h
    # collides with the right endpoint of the out block of k. that's because the
    # collision contributes to the uncertainty in the right endpoint of the in
    # block of the composition
    new_right += h.in_right - k.out_right
  end
  
  Exchanger(
    new_left,
    new_right,
    h.f_shift + k.f_shift,
    h.f_transit * k.f_transit,
    k.b_transit * h.b_transit
  )
end

function Base.show(io::IO, h::Exchanger)
  @printf(
    io,
    "[%0.3f, %0.3f] => [%0.3f, %0.3f]",
    map(mid, [h.in_left, h.in_right, h.out_left, h.out_right])...
  )
end

# === interval exchange cocycles

type Cocycle
  # we're using interval arithmetic to keep track of endpoints, so the intervals
  # of the interval exchange are called blocks to avoid confusion
  blocks::Array{Exchanger, 1}
end

# if you label the blocks 1, 2, 3... from left to right, apply the interval
# exchange, and read off the labels from left to right again, you get the list
# f_shuffle
function Cocycle{
  R <: Interval, S <: Any, T <: Integer
}(
  in_breaks::Array{R, 1}, f_transit::Array{S, 1}, f_shuffle::Array{T, 1}
)
  # find the lengths of the in blocks
  blocklengths = R[]
  last_b = 0
  for b in in_breaks
    push!(blocklengths, b - last_b)
    last_b = b
  end
  
  # add the break at zero
  pad_in_breaks = copy(in_breaks)
  unshift!(pad_in_breaks, 0)
  
  # find the breaks between the out blocks
  pad_out_breaks = cumsum([blocklengths[s] for s in f_shuffle])
  unshift!(pad_out_breaks, 0)
  
  # invert the permutation of the blocks
  b_shuffle = Array(T, length(f_shuffle))
  for (s, t) in enumerate(f_shuffle)
    b_shuffle[t] = s
  end
  
  # build the interval exchange, block by block
  blocks = Exchanger[]
  for (t, s) in enumerate(b_shuffle)
    push!(
      blocks,
      Exchanger(
        pad_in_breaks[t],
        pad_in_breaks[t+1],
        pad_out_breaks[s] - pad_in_breaks[t],
        f_transit[t],
        inv(f_transit[t])
      )
    )
  end
  
  Cocycle(blocks)
end

# compose an interval exchange cocyle with itself, roughly doubling the number
# of blocks (if the original interval exchange has n blocks, the new one will
# have 2n - 1)
function twostep(a::Cocycle)
  new_blocks = Exchanger[]
  for k in a.blocks
    # the number of comparisons here could be cut down significantly by taking
    # advantage of the fact that we know where the out block of k is, and
    # a.blocks is ordered by location of the in block
    for h in a.blocks
      piped = pipe(h, k)
      if piped != nothing
        push!(new_blocks, piped)
      end
    end
  end
  Cocycle(new_blocks)
end

function lamination(a::Cocycle, sym=nothing, depth=0)
  lam = Context[]
  
  # like the one in twostep, this loop could be made much more efficient
  for k in a.blocks
    for s in 1:(length(a.blocks) - 1)
      if !missed_connection(a.blocks[s], k) && !missed_connection(a.blocks[s+1], k)
        push!(
          lam,
          triangle_orbit(
            repeller(k.b_transit),
            repeller(a.blocks[s].f_transit),
            repeller(a.blocks[s+1].f_transit),
            50,
            sym,
            depth
          )
        )
        #=push!(
          lam,
          horoleaves(repeller(k.b_transit), repeller(a.blocks[s].f_transit), repeller(a.blocks[s+1].f_transit), 50)
        )=#
      end
    end
  end
  
  # ------
  
  out_order = sort(a.blocks, by=bl -> bl.out_left)
  for h in out_order
    for s in 1:(length(out_order) - 1)
      if !missed_connection(h, out_order[s]) && !missed_connection(h, out_order[s+1])
        push!(
          lam,
          triangle_orbit(
            repeller(h.f_transit),
            repeller(out_order[s+1].b_transit),
            repeller(out_order[s].b_transit),
            50,
            map(inv, sym),
            depth
          )
        )
      end
    end
  end
  
  compose(context(), lam...)
end

# === testing

include("regular.jl")

function lam_test()
  a = Cocycle(
    #[@interval(sin(0.3)), @interval(sin(0.6)), @interval(sin(1.2)), @interval(1)],
    [@interval(sin(0.1)), @interval(sin(0.5)), @interval(sin(1.4)), @interval(1)],
    generators(2),
    [4, 3, 2, 1]
  )
  
  for h in a.blocks
    println(h)
  end
  
  for i in 1:3
    a = twostep(a)
  end
  
  draw(
    PDF("lam_test.pdf", 10cm, 10cm),
    compose(context(), lamination(a, generators(2), 4), stroke("black"), linewidth(0.1))
  )
end

function arc_test()
  p = [1, 1im, -1im]
  q = [cis(0.2), cis(0.5), cis(1.3)]
  r = [cis(1.6), cis(3.4), cis(4.5)]
  draw(
    PDF("arc_test.pdf", 10cm, 10cm),
    compose(
      context(), linewidth(0.1),
      (
        context(0.1, 0.1, 0.8, 0.8),
        (
          context(), stroke("plum"),
          [horoleaves(p[t + 1], p[(t+1)%3 + 1], p[(t+2)%3 + 1], 50) for t in 0:2]...,
        ),
        (
          context(), stroke("lightsalmon"),
          [horoleaves(q[t + 1], q[(t+1)%3 + 1], q[(t+2)%3 + 1], 50) for t in 0:2]...,
        ),
        (
          context(), stroke("burlywood"),
          [horoleaves(r[t + 1], r[(t+1)%3 + 1], r[(t+2)%3 + 1], 50) for t in 0:2]...,
        ),
        (context(), circle(), fill("white"))
      ),
      (context(), rectangle(), fill("dimgrey"))
    )
  )
end

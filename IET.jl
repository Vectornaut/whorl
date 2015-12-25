module IET

using MPFI

# === utilities

# --- directions

@enum Direction forward=1 backward=-1

# --- interval arithmetic collision detection

type CollisionException <: Exception
  x::Interval
  y::Interval
end

function strict_lt(x::Interval, y::Interval)
  c = cmp(x, y)
  if c == 0
    throw(CollisionException(x, y))
  end
  c < 0
end

strict_gt(x::Interval, y::Interval) = !strict_lt(x, y)

# === interval exchange cocycles

# --- exchangers

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

# apply the exchanger to a point
function apply(h::Exchanger, x::Interval, dir::Direction=forward)
  if strict_lt(x, h.in_left) || strict_gt(x, h.in_right)
    return nothing
  end
  x + dir * h.f_shift
end

# compose the exchangers h and k, returning nothing if the out block of k
# doesn't overlap the in block of h
function pipe(h::Exchanger, k::Exchanger)
  # if out block of k doesn't overlap the in block of h, return nothing
  if strict_lt(h.in_right, k.out_left) || strict_gt(h.in_left, k.out_right)
    return nothing
  end
  
  # find the left endpoint of the in block of the composed exchanger
  new_left = k.in_left
  if strict_gt(h.in_left, k.out_left)
    new_left += h.in_left - k.out_left
  end
  
  # find the right endpoint of the in block of the composed exchanger
  new_right = k.in_right
  if strict_lt(h.in_right, k.out_right)
    new_right = k.in_right + h.in_right - k.out_right
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
    map(
      x -> convert(AbstractFloat, x),
      [h.in_left, h.in_right, h.out_left, h.out_right]
    )...
  )
end

# --- interval exchange cocycles

type Cocycle
  # we're using interval arithmetic to keep track of endpoints, so the intervals
  # of the interval exchange are called blocks to avoid confusion
  blocks::Array{Exchanger, 1}
end

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
end

# === testing

function test_routine()
  a = Cocycle(
    [map(sin, map(Interval, ["0.3", "0.6", "1.2"])); Interval("1")],
    [eye(2) for i in 1:4],
    [4, 3, 2, 1]
  )
  
  a2 = twostep(a)
  for h in a2.blocks
    println(h)
  end
end

end # module

include("poincare_disk.jl")

module IntervalExchange

using ValidatedNumerics, Compose, PoincaréDisk

export Cocycle, missed_connection, scancollect, twostep, lamination, foliage

# === exchangers

# an Exchanger is a piece of an interval exchange cocycle. it maps points in an
# "in" block to points in an "out" block by translation, discarding points that
# don't belong to the in block. it also applies a transition map to the
# local system sections above the block.
type Exchanger{R <: AbstractInterval}
  # the endpoints of the in and out blocks, specified using interval arithmetic.
  # the differences in_right - in_left and out_right - out_left will always be
  # equal, because the in block is mapped to the out block by translation.
  in_left::R
  in_right::R
  out_left::R
  out_right::R
  
  # the translation that sends the in block to the out block, always equal to
  # out_left - in_left
  f_shift::R
  
  # f_transit is the local system transition you pick up going forward through
  # the interval, and b_transit is the transition you pick up going backward.
  # f_transit and b_transit will always be inverses of one another.
  f_transit
  b_transit
  
  # the singularity at the right endpoint
  sing
  
  function Exchanger(in_left, in_right, f_shift, f_transit, b_transit, sing = nothing)
    new(
      in_left,            # in_left
      in_right,           # in_right
      in_left + f_shift,  # out_left
      in_right + f_shift, # out_right
      f_shift,            # f_shift
      f_transit,          # f_transit
      b_transit,          # b_transit
      sing                # sing
    )
  end
end

# interval arithmetic collision exception
type EndpointCollisionException <: Exception
  h::Exchanger
  k::Exchanger
end

# comparison function for sorting by in block
in_isless(h::Exchanger, k::Exchanger) =
  strictprecedes(h.in_left, k.in_left) && strictprecedes(h.in_right, k.in_right)

# comparison function for sorting by out block
out_isless(h::Exchanger, k::Exchanger) =
  strictprecedes(h.out_left, k.out_left) && strictprecedes(h.out_right, k.out_right)

# checks whether the out block of k hangs off to the right of the in block of h
function pre_hanging(h::Exchanger, k::Exchanger)
  # if we can't tell, throw an exception
  if !isdisjoint(h.in_right, k.out_right)
    throw(EndpointCollisionException(h, k))
  end
  
  strictprecedes(h.in_right, k.out_right)
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
function pipe{R <: AbstractInterval}(h::Exchanger{R}, k::Exchanger{R})
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
    
    # the right endpoint of the new block comes from to h
    new_sing = h.sing
  else
    # the right endpoint of the new block comes from k
    new_sing = k.sing
  end
  
  Exchanger{R}(
    new_left,
    new_right,
    h.f_shift + k.f_shift,
    h.f_transit * k.f_transit,
    k.b_transit * h.b_transit,
    new_sing
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

type Cocycle{R <: AbstractInterval}
  # we're using interval arithmetic to keep track of endpoints, so the intervals
  # of the interval exchange are called blocks to avoid confusion
  blocks_by_in::Array{Exchanger{R}, 1}
  blocks_by_out::Array{Exchanger{R}, 1}
end

# if you label the blocks 1, 2, 3... from left to right, apply the interval
# exchange, and read off the labels from left to right again, you get the list
# f_shuffle
function Cocycle{
  R <: AbstractInterval, S <: Any, T <: Integer
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
  blocks_by_in = Exchanger{R}[]
  for (t, s) in enumerate(b_shuffle)
    push!(
      blocks_by_in,
      Exchanger{R}(
        pad_in_breaks[t],
        pad_in_breaks[t+1],
        pad_out_breaks[s] - pad_in_breaks[t],
        f_transit[t],
        inv(f_transit[t])
      )
    )
  end
  
  # sort the block list by out block
  blocks_by_out = [blocks_by_in[s] for s in f_shuffle]
  
  Cocycle(blocks_by_in, blocks_by_out)
end

# scan through the in and out blocks of an interval exchange from left to right,
# applying functions to certain arrangements of blocks and collecting the
# results in order
# - when the in block of h connects to the out block of k, call thru_fn(h, k)
# - when the in blocks of left and right are adjacent, and both connect to the
#   out block of pivot, call f_fn(left, right, pivot)
# - when the out blocks of left and right are adjacent, and both connect to the
#   in block of pivot, call b_fn(left, right, pivot)
function scancollect(
  a::Cocycle,
  output_type::DataType;
  thru_fn::Union{Function, Void} = nothing,
  f_fn::Union{Function, Void} = nothing,
  b_fn::Union{Function, Void} = nothing
)
  output = output_type[]
  s = 1
  t = 1
  while true
    h = a.blocks_by_in[s]
    k = a.blocks_by_out[t]
    if thru_fn != nothing
      push!(output, thru_fn(h, k))
    end
    
    if s == length(a.blocks_by_in) && t == length(a.blocks_by_out)
      break
    elseif pre_hanging(h, k)
      s += 1
      if f_fn != nothing
        push!(output, f_fn(h, a.blocks_by_in[s], k))
      end
    else
      t += 1
      if b_fn != nothing
        push!(output, b_fn(k, a.blocks_by_out[t], h))
      end
    end
  end
  
  # return
  output
end

# compose an interval exchange cocyle with itself, roughly doubling the number
# of blocks (if the original interval exchange has n blocks, the new one will
# have 2n - 1)
function twostep{R <: AbstractInterval}(a::Cocycle{R})
  new_blocks = scancollect(a, Exchanger{R}, thru_fn = (h, k) -> pipe(h, k))
  new_blocks_by_in = sort(new_blocks, lt = in_isless)
  new_blocks_by_out = sort(new_blocks, lt = out_isless)
  Cocycle(new_blocks_by_in, new_blocks_by_out)
end

lamination{R <: AbstractInterval}(a::Cocycle{R}, sym, depth = 0) =
  scancollect(
    a, Context,
    thru_fn = (h, k) -> geodesic_orbit(
      repeller(k.b_transit),
      repeller(h.f_transit),
      sym, depth
    )
  )

foliage(a::Cocycle) =
  scancollect(
    a, Context,
    f_fn = (left, right, pivot) -> compose(
      context(),
      horotriangle(
        repeller(pivot.b_transit),
        repeller(left.f_transit),
        repeller(right.f_transit),
        69, 1/21, 4e-3
      ),
      stroke("coral"),
      linewidth(0.1mm)
    ),
    b_fn = (left, right, pivot) -> compose(
      context(),
      horotriangle(
        repeller(pivot.f_transit),
        repeller(right.b_transit),
        repeller(left.b_transit),
        69, 1/21, 4e-3
      ),
      stroke("deeppink"),
      linewidth(0.1mm)
    )
  )

end

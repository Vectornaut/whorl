include("poincare_disk.jl")

module IntervalExchange

using ValidatedNumerics, Compose, PoincaréDisk

export Cocycle, missed_connection, scancollect, twostep, Jump, FJump, BJump, abelianize

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
  
  # when we build an interval exchange cocycle from scratch, we label each
  # block with its in and out indices. by keeping track of this information
  # during composition, we can remember which blocks of the original cocycle the
  # composed exchanger starts and ends in
  orig_in::Integer
  orig_out::Integer
  
  # the singularity at the right endpoint
  sing
  
  function Exchanger(
    in_left,
    in_right,
    f_shift,
    f_transit,
    b_transit,
    orig_in,
    orig_out,
    sing = nothing
  )
    new(
      in_left,            # in_left
      in_right,           # in_right
      in_left + f_shift,  # out_left
      in_right + f_shift, # out_right
      f_shift,            # f_shift
      f_transit,          # f_transit
      b_transit,          # b_transit
      orig_in,            # orig_in
      orig_out,           # orig_out
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
    k.orig_in,
    h.orig_out,
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
        inv(f_transit[t]),
        t, s
      )
    )
  end
  
  # sort the block list by out block
  blocks_by_out = [blocks_by_in[s] for s in f_shuffle]
  
  # label the singularities
  sing = 1
  for start in 1:length(blocks_by_in)
    if blocks_by_in[start].sing == nothing
      # if we've found a new singularity, label all the blocks touching it
      s = start
      fwd = true
      while true
        # step to the next block around the singularity
        if fwd
          blocks_by_in[s].sing = sing
          s = b_shuffle[s] + 1
        else
          s = f_shuffle[s] - 1
        end
        
        # switch direction, unless we're going around an endpoint
        if s < 1
          s += 1
        elseif s > length(f_shuffle)
          s -= 1
        else
          fwd = !fwd
        end
        
        # if we're back where we started, break
        if (s, fwd) == (start, true)
          break
        end
      end
      
      # increment the label
      sing += 1
    end
  end
  
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
  thru_fn = nothing,
  f_fn = nothing,
  b_fn = nothing
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

# === abelianization

abstract Jump

Base.isless(p::Jump, q::Jump) = p.gap < q.gap

function jumpshear{T <: Number}(new_ln::Array{T}, old_ln::Array{T}, pivot::Array{T})
  new_sc = new_ln / det([new_ln pivot])
  old_sc = old_ln / det([old_ln pivot])
  [new_sc pivot] / [old_sc pivot]
end

type FJump <: Jump
  # the jump operator
  op
  
  # the singularity where the left and right blocks meet
  sing::Integer
  
  # the indices of the left, right, and pivot blocks in the original cocycle
  left::Integer
  right::Integer
  pivot::Integer
  
  # the stable lines
  left_stable
  right_stable
  pivot_stable
  
  # the square of the sine distance between the left and right stable lines
  gap
  
  function FJump{E <: Exchanger}(left::E, right::E, pivot::E)
    left_stable = stable(left.f_transit)
    right_stable = stable(right.f_transit)
    pivot_stable = stable(pivot.b_transit)
    new(
      jumpshear(left_stable, right_stable, pivot_stable), # op
      left.sing,                                          # sing
      left.orig_in,                                       # left
      right.orig_in,                                      # right
      pivot.orig_out,                                     # pivot
      left_stable,                                        # left_stable
      right_stable,                                       # right_stable
      pivot_stable,                                       # pivot_stable
      1 - abs(dot(left_stable, right_stable))             # gap
    )
  end
end

type BJump <: Jump
  # the jump operator
  op
  
  # the singularity where the left and right blocks meet
  sing::Integer
  
  # the indices of the left, right, and pivot blocks in the original cocycle
  left::Integer
  right::Integer
  pivot::Integer
  
  # the stable lines
  left_stable
  right_stable
  pivot_stable
  
  # the square of the sine distance between the left and right stable lines
  gap
  
  function BJump{E <: Exchanger}(left::E, right::E, pivot::E)
    left_stable = stable(left.b_transit)
    right_stable = stable(right.b_transit)
    pivot_stable = stable(pivot.f_transit)
    new(
      jumpshear(left_stable, right_stable, pivot_stable), # op
      left.sing,                                          # sing
      left.orig_out,                                      # left
      right.orig_out,                                     # right
      pivot.orig_in,                                      # pivot
      left_stable,                                        # left_stable
      right_stable,                                       # right_stable
      pivot_stable,                                       # pivot_stable
      1 - abs(dot(left_stable, right_stable))             # gap
    )
  end
end

# this function takes in two versions of the same SL(2,C) cocycle---orig, which
# was built from scratch using the Cocycle constructor, and iter, which has been
# iterated using twostep or a similar function. it returns the abelianized
# cocycle over the original interval exchange.
function abelianize(orig::Cocycle, iter::Cocycle)
  # build the abelianized cocycle, block by block
  blocks_by_in = typeof(first(orig.blocks_by_in))[]
  for bl in orig.blocks_by_in
    # we want to compute the abelianized holonomy around the following path:
    # - starting at the base segment, drive down the right lane of the left edge
    #   of bl's in block
    # - return to the base segment along the right lane of the left edge of bl's
    #   out block
    # - return to the starting point along the base segment
    # in comments, we'll refer to the right lane of the left edge of bl's in
    # block as y, the right lane of the left edge of bl's out block as x
    
    # check whether we'll be driving left or right when we return to the
    # starting point along the base segment
    leftward = strictprecedes(bl.in_left, bl.out_left)
    
    # set up a filter function that decides whether a given jump is on the way
    # from the x to y
    enroute = leftward ?
      # going left
      j -> if isa(j, FJump)
          bl.orig_in <= j.left && j.pivot < bl.orig_out
        elseif isa(j, BJump)
          bl.orig_in <= j.pivot && j.left < bl.orig_out
        end :
      # going right
      j -> if isa(j, FJump)
          bl.orig_in > j.left && j.pivot >= bl.orig_out
        elseif isa(j, BJump)
          bl.orig_in > j.pivot && j.left >= bl.orig_out
        end
    
    # compute the deviation from x to y by multiplying together all the
    # abelianization jumps on the way from x to y along the base interval. the
    # jumps come ordered from right to left, so we have to reverse the product
    # if the path from x to y goes to the right.
    ab_jumps = scancollect(iter, Jump, f_fn = FJump, b_fn = BJump)
    dev = prod([j.op for j in filter(enroute, ab_jumps)])
    if !leftward
      dev = inv(dev)
    end
    
    # compose the deviation with the parallel transport along the right lane of
    # the left edge of bl to get the abelianized holonomy
    ab_hol = dev * bl.f_transit
    
    # switch to the frame given by the forward- and backward-stable lines at y,
    # where the abelianized holonomy is diagonal
    stableframe = bl.orig_in == 1 ?
      begin
        # there's no jump over the left edge of the first in block, so the
        # forward- and backward-stable lines come from the leftmost in and out
        # blocks
        f_stable = stable(iter.blocks_by_in[1].f_transit)
        b_stable = stable(iter.blocks_by_out[1].b_transit)
        [f_stable b_stable]
      end :
      begin
        # find the jump over the left edge of the block and get the forward- and
        # backward lines from there
        edge = findfirst(j -> isa(j, FJump) && j.right == bl.orig_in, ab_jumps)
        f_stable = ab_jumps[edge].right_stable
        b_stable = ab_jumps[edge].pivot_stable
        stableframe = [f_stable b_stable]
      end
    diag_hol = stableframe \ ab_hol * stableframe
    
    # add the abelianized version of bl to the abelianized cocycle
    ab_bl = deepcopy(bl)
    ab_bl.f_transit = diag_hol
    ab_bl.b_transit = inv(diag_hol)
    push!(blocks_by_in, ab_bl)
  end
  
  # return
  Cocycle(
    blocks_by_in,
    sort(blocks_by_in, by = b -> b.orig_out)
  )
end

end

module Crawl

export findhome!, mapcollect, altcollect, tipcollect, CayleyCrawler, FreeCrawler, TileCrawler

const ROOT     =  3
const CLIMBER  =  0
const L_RUNNER = -1
const R_RUNNER =  1
const CORNER   =  2

abstract CayleyCrawler

function findhome!(crawler::CayleyCrawler, transit, base=eye(2))
  # find our way home
  if crawler.mode == ROOT
    crawler.home = base
  else
    crawler.home = base * transit[crawler.down+1]
  end
  
  # tell our shoots to find their ways home
  for sh in crawler.shoots
    findhome!(sh, transit, crawler.home)
  end
end

function mapcollect(f::Function, crawler::CayleyCrawler; prune = false)
  value = f(crawler.home)
  if prune && value == nothing
    collected = []
  else
    collected = [value]
  end
  if isempty(crawler.shoots)
    return collected
  else
    return vcat(collected, [mapcollect(f, sh, prune = prune) for sh in crawler.shoots]...)
  end
end

altcollect(f::Function, crawler::CayleyCrawler, include = true) =
  if isempty(crawler.shoots)
    return include ? [f(crawler.home)] : []
  else
    return vcat(
      include ? f(crawler.home) : [],
      [altcollect(f, sh, !include) for sh in crawler.shoots]...
    )
  end

function tipcollect(f::Function, crawler::CayleyCrawler; prune = false)
  if isempty(crawler.shoots)
    value = f(crawler.home)
    if prune && value == nothing
      return []
    else
      return [value]
    end
  else
    return vcat([tipcollect(f, sh, prune = prune) for sh in crawler.shoots]...)
  end
end

type FreeCrawler <: CayleyCrawler
  # the index of the generator represented by this edge of the Cayley graph,
  # oriented toward the root of the spanning tree (relevant in all modes but
  # root)
  down::Union{Integer, Void}
  
  # growth parameters
  mode::Integer   # specifies branching behavior
  ascent::Integer # how high to climb
  
  # the next edges out in the Cayley tree
  shoots::Array{FreeCrawler}
  
  # when the tile this edge ends on is centered, applying the home map centers
  # the root tile of the spanning tree
  home
  
  # a crawler for the free group on j generators, which has Schläfli symbol
  # {2j, ∞}
  FreeCrawler(j, ascent) = FreeCrawler(j, nothing, ROOT, ascent)
  
  function FreeCrawler(j, up, mode, ascent)
    # set parameters
    if mode != ROOT
      down = mod(up+j, 2j)
    else
      down = nothing
    end
    me = new(
      down,   # down
      mode,   # mode
      ascent, # ascent
      [],     # shoots
      nothing # home
    )
    
    # sprout new edges
    if mode == ROOT
      # the root sprouts climbers in all directions
      for m in 0:(2j - 1)
        push!(me.shoots, FreeCrawler(j, m, CLIMBER, ascent))
      end
    elseif mode == CLIMBER && ascent > 1
      # a climber sprouts climbers in every direction but down
      for m in 1:(2j - 1)
        push!(me.shoots, FreeCrawler(j, down+m, CLIMBER, ascent-1))
      end
    end
    
    # return
    me
  end
end

type TileCrawler <: CayleyCrawler
  # the index of the generator represented by this edge of the Cayley graph,
  # oriented toward the root of the spanning tree (relevant in all modes but
  # root)
  down::Union{Integer, Void}
  
  # growth parameters
  mode::Integer               # specifies branching behavior
  range::Union{Integer, Void} # how far to run (relevant in the runner modes)
  ascent::Integer             # how high to climb
  
  # the next edges out in the spanning tree
  shoots::Array{TileCrawler}
  
  # when the tile this edge ends on is centered, applying the home map centers
  # the root tile of the spanning tree
  home
  
  # a crawler for the group with Schläfli symbol {2j, 2k}
  TileCrawler(j, k, ascent) = TileCrawler(j, k, nothing, ROOT, nothing, ascent)
  
  function TileCrawler(j, k, up, mode, range, ascent)
    # set parameters
    if mode != ROOT
      down = mod(up+j, 2j)
    else
      down = nothing
    end
    me = new(
      down,   # down
      mode,   # mode
      range,  # range
      ascent, # ascent
      [],     # shoots
      nothing # home
    )
    
    # sprout new edges
    if mode == ROOT
      # the root sprouts climbers in all directions
      for m in 0:(2j - 1)
        push!(me.shoots, TileCrawler(j, k, m, CLIMBER, nothing, ascent))
      end
    elseif mode == CLIMBER
      # a climber sprouts a left-runner in the left position, a right-runner in
      # the right position, and more climbers in between
      push!(me.shoots, TileCrawler(j, k, down-1, L_RUNNER, k-2, ascent))
      push!(me.shoots, TileCrawler(j, k, down+1, R_RUNNER, k-2, ascent))
      if ascent > 1
        for m in 2:(2j - 2)
          push!(me.shoots, TileCrawler(j, k, down+m, CLIMBER, nothing, ascent-1))
        end
      end
    elseif mode == L_RUNNER
      # a left-runner's behavior depends on its range
      if range > 1
        # if it's not at the end of its range, it sprouts another left-runner in
        # the left position
        push!(me.shoots, TileCrawler(j, k, down-1, L_RUNNER, range-1, ascent))
      else
        # if it's at the end of its range, it sprouts a corner in the left position
        push!(me.shoots, TileCrawler(j, k, down-1, CORNER, nothing, ascent))
      end
      
      # it sprouts a right-runner in the right position, and climbers in between
      if ascent > 1
        push!(me.shoots, TileCrawler(j, k, down+1, R_RUNNER, k-2, ascent-1))
        for m in 2:(2j - 2)
          push!(me.shoots, TileCrawler(j, k, down+m, CLIMBER, nothing, ascent-1))
        end
      end
    elseif mode == R_RUNNER
      # a right-runner's behavior depends on its range
      if range > 1
        # if it's not at the end of its range, it sprouts another right-runner
        # in the right position
        push!(me.shoots, TileCrawler(j, k, down+1, R_RUNNER, range-1, ascent))
      end
      # if it's at the end of its range, it sprouts nothing in the right position
      
      # it sprouts a left-runner in the left position, and climbers in between
      if ascent > 1
        push!(me.shoots, TileCrawler(j, k, down-1, L_RUNNER, k-2, ascent-1))
        for m in 2:(2*j - 2)
          push!(me.shoots, TileCrawler(j, k, down+m, CLIMBER, nothing, ascent-1))
        end
      end
    elseif mode == CORNER
      # a corner sprouts nothing in the left position, a left-runner in the
      # second-left position, a right-runner in the right position, and climbers
      # in between
      if ascent > 1
        push!(me.shoots, TileCrawler(j, k, down-2, L_RUNNER, k-2, ascent-1))
        push!(me.shoots, TileCrawler(j, k, down+1, R_RUNNER, k-2, ascent-1))
        for m in 2:(2j - 3)
          push!(me.shoots, TileCrawler(j, k, down+m, CLIMBER, nothing, ascent-1))
        end
      end
    end
    
    # return
    me
  end
end

end

module Lamination

export
  SOLID,
  HORO,
  LaminationTheme,
  tacos,
  bone,
  shell,
  polygon_penciller,
  polygon_inker,
  triangle_inker,
  IdealPolygon,
  orbiter,
  render

using
  Colors,
  Compose,
  PoincaréDisk,
  Crawl,
  ValidatedNumerics,
  IntervalExchange,
  Caterpillar,
  Square,
  PunkdTorus

# === themes

# fill types
const SOLID = 0
const HORO = 1

type LaminationTheme
  leafcolor  # the color of the leaves
  fillcolor  # the colors of the complementary triangles
  fillstyle  # solid or horocyclic
  checkcolor # the color of the fundamental domain checkers
  diskcolor  # the color of the background disk
end

const tacos = LaminationTheme(
  nothing, # leafcolor
  [
    RGB(255/255, 1/255, 73/255),
    RGB(255/255, 121/255, 1/255),
    RGB(255/255, 210/255, 0/255),
    RGB(0/255, 200/255, 146/255)
  ],       # fillcolor
  SOLID,   # fillstyle
  nothing, # checkcolor
  nothing  # diskcolor
)

const bone = LaminationTheme(
  RGB(172/255, 146/255, 122/255),       # leafcolor
  fill(RGB(99/255, 96/255, 83/255), 4), # fillcolor
  HORO,                                 # fillstyle
  nothing,                              # checkcolor
  RGB(39/255, 36/255, 33/255)           # diskcolor
)

const shell = LaminationTheme(
  RGB(0/255, 30/255, 140/255),           # leafcolor
  fill(RGB(0/255, 150/255, 173/255), 4), # fillcolor
  HORO,                                  # fillstyle
  RGB(235/255, 233/255, 229/255),        # checkcolor
  RGB(1, 1, 1)                           # diskcolor
)

polygon_penciller(theme::LaminationTheme) =
  (verts, sing) -> compose(
    context(),
    ideal_edges(verts...),
    stroke(theme.leafcolor)
  )

polygon_inker(theme::LaminationTheme) =
  (verts, sing) -> compose(
    context(),
    ideal_path(verts...),
    fill(theme.checkcolor)
  )

function triangle_inker(theme::LaminationTheme)
  if theme.fillstyle == SOLID
    return (verts, sing) -> compose(
      context(),
      ideal_path(verts...),
      fill(theme.fillcolor[sing])
    )
  elseif theme.fillstyle == HORO
    return (verts, sing) -> compose(
      context(),
      horotriangle(verts..., 69, 1/21, 4e-3),
      stroke(theme.fillcolor[sing])
    )
  end
end

# === ideal polygons

type IdealPolygon
  sing::Integer
  verts
end

function startblock(orig_out, a::Cocycle)
  index = findfirst(bl -> bl.orig_out == orig_out, a.blocks_by_out)
  a.blocks_by_out[index]
end

function endblock(orig_out, a::Cocycle)
  index = findlast(bl -> bl.orig_out == orig_out, a.blocks_by_out)
  a.blocks_by_out[index]
end

# find the triangle formed by the forward- and backward-stable lines at a jump
# in a cocycle
triangulate(j::Jump) =
  IdealPolygon(
    j.sing,
    [planeproj(v) for v in [j.pivot_stable, j.right_stable, j.left_stable]]
  )

# foliate the translation surface carrying `loc` at the given angle, pull the
# leaves tight with respect to the hyperbolic structure described by `loc`, and
# return the complementary triangles of the resulting geodesic lamination. we
# approximate the lamination by doubling the first return cocycle `depth` times.

function triangulate{R <: AbstractInterval}(angle::R, loc::CaterpillarLocSys, depth::Integer; verbose = false)
  # build and evolve cocycle
  orig = Caterpillar.cocycle(angle, loc)
  iter = power_twostep(orig, depth, verbose = verbose)
  
  # find the widest triangle for each singularity
  b_jumps = scancollect(iter, Jump, b_fn = BJump)
  widest = Jump[]
  for sing in 1:4
    push!(widest, maximum(filter(j -> j.sing == sing, b_jumps)))
  end
  
  # return
  [triangulate(j) for j in widest]
end

function triangulate{R <: AbstractInterval}(angle::R, loc::PunkdTorusLocSys, depth::Integer; verbose = false)
  # set up cocycle
  orig = PunkdTorus.cocycle(angle, loc)
  iter = power_twostep(orig, depth, verbose = verbose)
  
  # find the triangle to the right of the critical leaf rising out of the
  # northwest puncture
  nw_verts = [planeproj(v) for v in [
    loc.punks[1],
    orig.blocks_by_out[2].b_transit * stable(startblock(2, iter).b_transit),
    stable(first(iter.blocks_by_in).f_transit)
  ]]
  
  # find the triangle to the left of the critical leaf rising out of the
  # southeast puncture
  se_verts = [planeproj(v) for v in [
    loc.punks[3],
    stable(last(iter.blocks_by_in).f_transit),
    orig.blocks_by_out[1].b_transit * stable(endblock(1, iter).b_transit)
  ]]
  
  # return
  [IdealPolygon(1, nw_verts), IdealPolygon(2, se_verts)]
end

# given an ideal polygon, return the function that takes a möbius transformation
# m, applies it to the polygon, and draws the result with the given drawing
# function
orbiter(p::IdealPolygon, eps, draw, shift = eye(2); diam = [2, 3]) =
  m -> begin
    verts = [möbius_map(shift*m, v) for v in p.verts]
    if eps == nothing || abs2(verts[diam[1]] - verts[diam[2]]) > eps*eps
      return draw(verts, p.sing)
    else
      return compose(context())
    end
  end

# === rendering

# foliate the translation surface carrying `loc` at the given angle, pull the
# leaves tight with respect to the hyperbolic structure described by `loc`, and
# draw the resulting geodesic lamination. we approximate the lamination by
# finding its complementary triangles and tiling them using the given crawler
# for the holonomy group. if a frame number is specified, render an
# appropriately named bitmap to be used as a frame of a movie. otherwise, render
# a PDF or SVG test frame.
function render{R <: AbstractInterval}(
  angle::R,
  loc,
  crawler,
  eps,
  theme;
  frame = nothing,
  center = nothing,
  svg = false,
  verbose = false
)
  # get complementary triangles
  triangles = triangulate(angle, loc, 4, verbose = verbose)
  
  # do centering, if requested
  if center == nothing
    shift = eye(2)
  else
    shift = pts_to_pts(
      triangles[center].verts...,
      cis(3π/6), cis(7π/6), cis(11π/6)
    )
  end
  
  # start a list of layers
  layers = []
  
  # draw leaves
  if theme.leafcolor != nothing
    clip = compose(context(), circle(), stroke("white"), linewidth(0.1mm), fill(nothing))
    push!(layers, clip)
    
    leaves = vcat([
      mapcollect(orbiter(t, eps, polygon_penciller(theme), shift), crawler)
      for t in triangles
    ]...)
    leaf_gp = compose(context(), leaves..., linewidth(0.1mm), fill(nothing))
    push!(layers, leaf_gp)
  end
  
  # fill complementary triangles
  if theme.fillcolor != nothing
    fills = vcat([
      mapcollect(orbiter(t, eps, triangle_inker(theme), shift), crawler)
      for t in triangles
    ]...)
    if theme.fillstyle == SOLID
      fill_gp = compose(context(), fills...)
    elseif theme.fillstyle == HORO
      fill_gp = compose(context(), fills..., linewidth(0.1mm), fill(nothing))
    end
    push!(layers, fill_gp)
  end
  
  # draw fundamental domain checkers
  if theme.checkcolor != nothing && isa(loc, PunkdTorusLocSys)
    fund = IdealPolygon(1, [planeproj(v) for v in loc.punks])
    checks = altcollect(orbiter(fund, eps, polygon_inker(theme), shift), crawler)
    check_gp = compose(context(), checks...)
    push!(layers, check_gp)
  end
  
  # draw background
  if theme.diskcolor != nothing
    disk = compose(context(), circle(), fill(theme.diskcolor), stroke(nothing))
    push!(layers, disk)
  end
  
  # render
  picture = compose(context(), layers...)
  if frame == nothing
    if svg
      draw(SVG("triangle_test.svg", 7cm, 7cm), picture)
    else
      draw(PDF("triangle_test.pdf", 7cm, 7cm), picture)
    end
  else
    draw(PNG(@sprintf("triangle_mov/frame%02i.png", frame), 500px, 500px), picture)
  end
end

end
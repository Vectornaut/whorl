module PunkdTorus

using
  ValidatedNumerics,
  IntervalExchange,
  Square,
  Crawl,
  PoincaréDisk,
  Compose
  ##Examples

export PunkdTorusLocSys

type PunkdTorusLocSys <: SquareLocSys
  n_transit
  e_transit
  s_transit
  w_transit
  
  function PunkdTorusLocSys(length, twist)
    # the mobius transformation
    # -i --> -1
    #  1 -->  0
    #  i -->  1
    flatten = [[-im, 1] [im, 1]]
    
    # see Parker and Parkkonen's "Coordinates for Quasi-Fuchsian Punctured Torus
    # Space", Figure 1.2. our `arcleft` is the inverse of their S', and our
    # `expand` is their T.
    arcleft = [
      -cosh(length)        cosh(length) - 1;
       cosh(length) + 1   -cosh(length)
    ]
    expand = [
       cosh(twist/2)*coth(length/2)   -sinh(twist/2);
      -sinh(twist/2)                   cosh(twist/2)*tanh(length/2)
    ]
    
    # return
    new(
      inv(flatten) * arcleft * flatten,
      inv(flatten) * expand * flatten,
      inv(flatten) * inv(arcleft) * flatten,
      inv(flatten) * inv(expand) * flatten
    )
  end
end

function punkd_torus_loc_sys(length, twist)
  # the mobius transformation
  # -i --> -1
  #  1 -->  0
  #  i -->  1
  flatten = [[-im, 1] [im, 1]]
  
  # see Parker and Parkkonen's "Coordinates for Quasi-Fuchsian Punctured Torus
  # Space", Figure 1.2. our `arcleft` is the inverse of their S', and our
  # `expand` is their T.
  arcleft = [
    -cosh(length)        cosh(length) - 1;
     cosh(length) + 1   -cosh(length)
  ]
  expand = [
     cosh(twist/2)*coth(length/2)   -sinh(twist/2);
    -sinh(twist/2)                   cosh(twist/2)*tanh(length/2)
  ]
  
  # return
  RectangleLocSys{AbstractInterval}(
    @interval(1), @interval(1),
    inv(flatten) * arcleft * flatten,
    inv(flatten) * expand * flatten
  )
end

type PunkJump <: Jump
  # the singularity where the left and right blocks meet
  sing::Integer
  
  # the stable lines
  left_stable
  right_stable
  pivot_stable
  
  function PunkJump{E <: IntervalExchange.Exchanger}(sing, left::E, cw_transit)
    left_stable = stable(left.f_transit)
    new(
      ##left.sing,                 # sing
      sing,
      left_stable,               # left_stable
      cw_transit * left_stable,  # right_stable
      eigvecs(cw_transit)[:,1],  # pivot_stable
    )
  end
  
  PunkJump(sing, left_stable, right_stable, pivot_stable) =
    new(sing, left_stable, right_stable, pivot_stable)
end

function render{R <: AbstractInterval}(
  angle::R,
  loc,
  crawler,
  eps,
  theme;
  frame = nothing,
  ##center = nothing,
  svg = false
)
  # set up cocycle
  orig = cocycle(angle, loc)
  
  # evolve cocycle
  iter = orig
  for i in 1:6
    print("  Step $i\n  ")
    iter = @time(twostep(iter))
    println("    $(length(iter.blocks_by_in)) blocks")
  end
  
  # write down a wide, convenient triangle for each singularity
  left_jump = PunkJump(
    1,
    orig.blocks_by_out[1].b_transit * stable(iter.blocks_by_out[findlast(bl -> bl.orig_out == 1, iter.blocks_by_out)].b_transit),
    stable(last(iter.blocks_by_in).f_transit),
    eigvecs(loc.e_transit * loc.n_transit * loc.w_transit * loc.s_transit)[:,1]
  )
  right_jump = PunkJump(
    2,
    stable(first(iter.blocks_by_in).f_transit),
    orig.blocks_by_out[2].b_transit * stable(iter.blocks_by_out[findfirst(bl -> bl.orig_out == 2, iter.blocks_by_out)].b_transit),
    eigvecs(loc.s_transit * loc.e_transit * loc.n_transit * loc.w_transit)[:,1]
  )
  corner1 = PunkJump(
    3,
    eigvecs(loc.e_transit * loc.n_transit * loc.w_transit * loc.s_transit)[:,1],
    eigvecs(loc.s_transit * loc.e_transit * loc.n_transit * loc.w_transit)[:,1],
    eigvecs(loc.w_transit * loc.s_transit * loc.e_transit * loc.n_transit)[:,1]
  )
  corner2 = PunkJump(
    4,
    eigvecs(loc.w_transit * loc.s_transit * loc.e_transit * loc.n_transit)[:,1],
    eigvecs(loc.n_transit * loc.w_transit * loc.s_transit * loc.e_transit)[:,1],
    eigvecs(loc.e_transit * loc.n_transit * loc.w_transit * loc.s_transit)[:,1]
  )
  jumps = [left_jump, right_jump, corner1, corner2]
  
  # do centering, if requested
  shift = eye(2) ## FILL IN LATER
  
  # start a list of layers
  layers = []
  
  # everything below this line is pretty much verbatim from the render function
  # in examples.jl
  # ---------
  
  # draw leaves
  if theme.leafcolor != nothing
    clip = compose(context(), circle(), stroke("white"), linewidth(0.1mm), fill(nothing))
    push!(layers, clip)
    
    leaves = vcat([mapcollect(Examples.leaf_orbiter(j, eps, theme, shift), crawler) for j in jumps]...)
    leaf_gp = compose(context(), leaves..., linewidth(0.1mm), fill(nothing))
    push!(layers, leaf_gp)
  end
  
  # fill complementary triangles
  if theme.fillcolor != nothing
    triangles = vcat([mapcollect(Examples.fill_orbiter(j, eps, theme, shift), crawler) for j in jumps]...)
    if theme.fillstyle == Examples.SOLID
      triangle_gp = compose(context(), triangles...)
    elseif theme.fillstyle == Examples.HORO
      triangle_gp = compose(context(), triangles..., linewidth(0.1mm), fill(nothing))
    end
    push!(layers, triangle_gp)
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
      draw(SVG("punkd_test.svg", 7cm, 7cm), picture)
    else
      draw(PDF("punkd_test.pdf", 7cm, 7cm), picture)
    end
  else
    draw(PNG(@sprintf("punkd_mov/frame%02i.png", frame), 500px, 500px), picture)
  end
end

function lam_test(angle = @interval(0.5))
  loc = punkd_torus_loc_sys(1,0)
  crawler = FreeCrawler(2,5);
  findhome!(crawler, [loc.n_transit, loc.e_transit, loc.s_transit, loc.w_transit])
  render(angle, loc, crawler, 1e-3, Examples.tacos)
end

end

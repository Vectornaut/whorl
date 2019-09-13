include("cayley_crawler.jl")
include("regular.jl")
include("interval_exchange.jl")
include("caterpillar.jl")
include("square.jl")
include("punkd_torus.jl")
include("quasicrystal.jl")
include("lamination.jl")

module Examples

using
  Printf,
  Gadfly,
  DataFrames,
  Colors,
  Compose,
  Main.PoincaréDisk,
  Main.Crawl,
  LinearAlgebra,
  ValidatedNumerics,
  Main.IntervalExchange,
  Main.Caterpillar,
  Main.Square,
  Main.PunkdTorus,
  Main.Quasicrystal,
  Main.Lamination

import Cairo, Fontconfig, Main.Regular

# === basic drawing

function farey_tiles(; eps = 1e-3, theme = shell, svg = false)
  # write down symmetry group generators in the upper half-plane model
  farey_generators = [[1 2; 0 1], [1 0; 2 1], [1 -2; 0 1], [1 0; -2 1]]
  
  # write down the mobius transformation
  # -i --> -1
  #  1 -->  0
  #  i -->  1
  flatten = [[-im, 1] [im, 1]]
  
  # set up crawler
  crawler = FreeCrawler(2, 5)
  findhome!(crawler, [inv(flatten) * g * flatten for g in farey_generators])
  
  # start a list of layers
  layers = []
  
  # draw triangle lifts
  base = IdealPolygon(1, [1, im, -1])
  leaves = mapcollect(orbiter(base, eps, polygon_penciller(theme)), crawler)
  leaf_layer = compose(context(), leaves..., stroke(theme.leafcolor), linewidth(0.1mm), fill(nothing))
  fills = mapcollect(orbiter(base, eps, polygon_inker(theme)), crawler)
  fill_layer = compose(context(), fills..., fill(theme.checkcolor), stroke(nothing))
  push!(layers, leaf_layer, fill_layer)
  
  # draw background
  if theme.diskcolor != nothing
    disk = compose(context(), Compose.circle(), fill(theme.diskcolor), stroke(nothing))
    push!(layers, disk)
  end
  
  # render
  picture = compose(context(), layers...)
  if svg
    picture |> SVG("farey_tiles.svg", 7cm, 7cm)
  else
    picture |> PDF("farey_tiles.pdf", 7cm, 7cm)
  end
end

# === abelianization

# print a nicely formatted complex number
function prettyprint(z::Number)
  @printf("%8.3f %c %7.3fim", real(z), imag(z) < 0 ? '-' : '+', abs(imag(z)))
end

# print a nicely formatted complex matrix
function prettyprint(m::Matrix, indent = 0)
  for i in 1:size(m, 1)
    print(" "^indent)
    for j in 1:size(m, 2)
      prettyprint(m[i,j])
      print("   ")
    end
    println()
  end
  println()
end

# given a matrix cocycle over an interval exchange, print the block translations
# of the interval exchange along with the transition maps of the cocycle
function printcocycle(a::Cocycle)
  i = 0
  for bl in a.blocks_by_in
    i += 1
    @printf("Block %-5d", i)
    println("$bl\n")
    prettyprint(bl.f_transit, 4)
  end
end

# print an "almost flat caterpillar" cocycle and its abelianization
function abelianization_ex()
  # build and evolve cocycle
  loc = almost_flat_caterpillar(Regular.generators(4, 4))
  orig = Caterpillar.cocycle(@interval(3π/4 + 1//11), loc)
  iter = power_twostep(orig, 4)
  
  # abelianize
  ab = IntervalExchange.abelianize(orig, iter)
  
  # output
  println("=== original cocyle\n")
  printcocycle(orig)
  println("=== abelianized cocyle\n")
  printcocycle(ab)
end

# === shear parameter plots

# compute the shear parameters of an SL(2,C) cocycle
function shears(orig)
  # evolve cocycle
  iter = power_twostep(orig, 4)
  
  # abelianize
  ab = IntervalExchange.abelianize(orig, iter)
  [bl.f_transit[1,1] for bl in ab.blocks_by_in]
end

# find the Masur polygon of a deflated genus-2 surface
function render_deflation(
  angle::R,
  g;
  theme = tacos,
  half = false
) where R <: AbstractInterval
  # build cocycle and compute shears along first-return paths
  loc = almost_flat_caterpillar(g)
  orig = Caterpillar.cocycle(angle, loc)
  x = shears(orig)
  
  # compute shears along odd paths around polygon sides
  y = fill(NaN + NaN*im, 13)
  y[4] = inv(x[3])*x[2]*x[7]*inv(x[6])*x[5]
  y[5] = inv(x[4])*x[3]*inv(x[1])*inv(x[13])*x[6]
  y[6] = inv(x[5])*x[4]*inv(x[2])*x[1]*x[7]
  y[7] = inv(x[6])*x[5]*inv(x[4])*inv(x[9])*x[8]
  y[8] = inv(x[7])*x[12]*inv(x[10])*x[9]
  y[9] = inv(x[8])*x[13]*inv(x[11])*x[10]
  y[10] = inv(x[9])*x[8]*x[13]*inv(x[12])*x[11]
  y[13] = inv(x[12])*x[11]*inv(x[10])*inv(x[3])*x[2]*inv(x[1])
  y[1] = y[4]
  y[2] = y[5]
  y[3] = y[6]
  y[11] = y[8]
  y[12] = y[9]
  
  # find side runs of Masur polygon
  runs = hcat(
    [mid(block.in_right - block.in_left) for block in orig.blocks_by_in],
    [-log(-real(shear)) for shear in y]
  )
  
  # list vertices
  f_shuffle = [7, 12, 11, 10, 9, 8, 13, 6, 5, 4, 3, 2, 1]
  runsaround = [
    [runs[t,:] for t in 1:length(y)];
    [-runs[t,:] for t in reverse(f_shuffle[2:end])]
  ]
  vtc = map(p -> tuple(p...), cumsum(runsaround))
  pushfirst!(vtc, (0, 0))
  
  # compute dimensions
  width = mid(orig.blocks_by_in[end].in_right)
  bottom = min([v[2] for v in vtc]...)
  height = max([v[2] for v in vtc]...) - bottom
  
  # return
  if half
    desired_vtc = vcat(vtc[4:10], vtc[17:23])
  else
    desired_vtc = vtc
  end
  compose(
    context(units = UnitBox(0, bottom, width, height)),
    (context(), polygon(desired_vtc)),
    candystripes(angle, loc, 3, theme)
  )
end

function deflation_movie(; ascent = 2, eps = 1e-3, theme = tacos, testframe = true, center = nothing, svg = false, verbose = true)
  angle = @interval(3π/4 + 1//11)
  crawler = TileCrawler(4, 4, ascent)
  if testframe
    n = 1
    timepoints = [0]
  else
    n = 25
    timepoints = 1:n
  end
  for t in timepoints
    println(testframe ? "Test frame" : "Frame $t")
    u = easing(t//n)
    
    # find Dehn-twisted holonomies
    transit = Regular.generators(4, 4)
    frame = [[1, 1] [-1, 1]] * [[im, 0] [0, -im]]
    twist = frame * Diagonal([1/exp(u/2), exp(u/2)]) * inv(frame)
    for s in 1:3
      transit[s] = twist*transit[s]
    end
    
    # set up crawler
    findhome!(crawler, [transit; [inv(t) for t in transit]])
    
    # draw frame
    lam_pic = @time(Lamination.render(angle, almost_flat_caterpillar(transit), crawler, eps, theme, center = center, verbose = verbose))
    defl_pic = @time(render_deflation(angle, transit, theme = theme))
    frame = compose(context(),
      compose(context(0.05, 0.025, 0.9, 0.45), lam_pic),
      compose(context(0.05, 0.525, 0.9, 0.45), defl_pic)
    )
    if testframe
      if svg
        draw(SVG("deflation_test.svg", 7cm, 14cm), frame)
      else
        draw(PDF("deflation_test.pdf", 7cm, 14cm), frame)
      end
    else
      draw(PNG(@sprintf("deflation_mov/frame%02i.png", t), 500px, 1000px), frame)
    end
  end
end

# to compose the final image from pre-poster.svg:
# 1) reflect the universal cover and the candy-striped Masur polygon vertically
# 2) reflect the Masur polygon and its candy stripes horizontally
# 3) separate the four colors into layers
# 4) clip each separated candy stripe layer with a copy of the Masur polygon
# pantone color assignments:
#   strawberry = Hexachrome Magenta U
#   orange = 1505 U
#   gold = 108 U
#   aqua = 3265 U
function bridges_poster()
  # set parameters. printed size of universal cover is about 812 points, and Kid
  # Icarus can print lines as thin as 1 point reliably, so eps should be below
  # 1.23e-3
  angle = @interval(3π/4 + 1//11)
  ascent = 4
  eps = 4e-4
  
  # set up crawler
  transit = Regular.generators(4, 4)
  crawler = TileCrawler(4, 4, ascent)
  findhome!(crawler, [transit; [inv(t) for t in transit]])
  
  # draw poster
  lam_pic = @time(Lamination.render(angle, almost_flat_caterpillar(transit), crawler, eps, tacos, verbose = true))
  defl_pic = @time(render_deflation(angle, transit, theme = tacos, half = true))
  full_w = 558.8;  full_h = 838.2
  lam_x = 82.550;  lam_y = 352.350; lam_side = 393.700
  defl_x = 36.001; defl_y = 39.430; defl_w = 486.797;  defl_h = 232.139
  poster = compose(context(),
    compose(context(lam_x/full_w, 1 - (lam_y + lam_side)/full_h, lam_side/full_w, lam_side/full_h), lam_pic),
    compose(context(defl_x/full_w, 1 - (defl_y + defl_h)/full_h, defl_w/full_w, defl_h/full_h), defl_pic)
  )
  draw(SVG("pre-poster.svg", 16inch, 24inch), poster)
end

# `range` doesn't work with Interval objects, so here's a slapdash replacement
function grid(start, fin, res::Integer)
  map(u -> (1-u)*start + u*fin, [@interval(t//(res-1)) for t in 0:res-1])
end

# compute shear parameters of a family of cocycles over a range of angles, with
# the given resolution. the argument `cyc` should be a function that takes an
# angle, of type AbstractInterval, and returns a cocycle. the range of angles
# goes roughly from π/2 to π. we use rational approximations for π/2 and π that
# are only accurate to about one part in 114, making it unlikely that we'll get
# within machine precision of a saddle connection at any reasonable resolution.
function shear_data(cyc, res::Integer)
  angles = grid(@interval(358//114), @interval(180//114), res)
  angle_col = []
  block_col = []
  real_shear_col = []
  imag_shear_col = []
  
  bunchsize = max(div(length(angles), 10), 10)
  i = 0
  while i < length(angles)
    println("Points $(i + 1) through $(min(i + bunchsize, length(angles)))")
    @time(
      for _ in 1:bunchsize
        i += 1
        if i > length(angles)
          break
        end
        
        x = shears(cyc(angles[i]))
        append!(angle_col, fill(mid(angles[i]), length(x)))
        append!(block_col, collect(1:length(x)))
        append!(real_shear_col, map(real, x))
        append!(imag_shear_col, map(imag, x))
      end
    )
  end
  
  DataFrame(
    Any[angle_col, block_col, real_shear_col, imag_shear_col],
    map(Symbol, ["angle", "block", "real(shear)", "imag(shear)"])
  )
end

function shear_plot(data)
  # scale for angles
  scale = Scale.x_continuous(minvalue = Float64(π/2), maxvalue = Float64(π))
  
  # ticks at the shortest saddle connections
  saddle_ticks = Guide.xticks(
    ticks = [
      map(
        t -> atan(t) + π/2,
        [0, 1/4, 1/3, 1/2, 2/3, 3/4, 1, 4/3, 3/2, 2/1, 3/1, 4/1]
      )
      π
    ],
    label = false
  )
  
  # themes
  real_theme = Theme(
    default_color = tacos.fillcolor[1],
    point_size = 0.5mm,
    highlight_width = 0mm,
  )
  imag_theme = Theme(
    default_color = tacos.fillcolor[4],
    point_size = 0.5mm,
    highlight_width = 0mm,
  )
  
  # plot
  real_layer = layer(x = "angle", y = "real(shear)", Geom.point, real_theme)
  imag_layer = layer(x = "angle", y = "imag(shear)", Geom.point, imag_theme)
  plot(
    data,
    ygroup = "block",
    Geom.subplot_grid(
      real_layer, imag_layer,
      scale, saddle_ticks,
      free_y_axis = true
    ),
    Guide.xlabel("angle"),
    Guide.ylabel("shear <b><i>by</i></b> block")
  )
end

traceless(h::Number, x::Number, y::Number) = [h x; y -h]

expconj(t) =
  tup -> begin
    g, a = tup
    exp(t*a)*g*exp(-t*a)
  end

function shear_ex(cyc, transit, perturbation; highres = false)
  # no perturbation
  println("=== no perturbation\n")
  cyc_no = cyc(transit)
  data_no = shear_data(cyc_no, highres ? 300 : 18)
  p_no = shear_plot(data_no)
  page_height = data_no.block[end]*3cm + 1cm # data_no.block[end] gives the number of blocks in the cocycle
  draw(PDF("no-perturbation.pdf", 30cm, page_height), p_no)
  println()
  
  # small perturbation
  println("=== small perturbation\n")
  transit_sm = map(expconj(0.01), zip(transit, perturbation))
  data_sm = shear_data(cyc(transit_sm), highres ? 300 : 18)
  p_sm = shear_plot(data_sm)
  draw(PDF("small-perturbation.pdf", 30cm, page_height), p_sm)
  println()
  
  # large perturbation
  println("=== large perturbation\n")
  transit_lg = map(expconj(0.1), zip(transit, perturbation))
  data_lg = shear_data(cyc(transit_lg), highres ? 300 : 18)
  p_lg = shear_plot(data_lg)
  draw(PDF("large-perturbation.pdf", 30cm, page_height), p_lg)
end

function square_shear_ex(k; highres = false)
  cyc = transit -> begin
    loc = GenSquareLocSys(transit...)
    angle -> Square.cocycle(angle, loc)
  end
  perturbation = Matrix[
    traceless(10, 0, 0),
    traceless(0, 10, 0)
  ]
  shear_ex(cyc, Regular.generators(2, k), perturbation, highres = highres)
end

function caterpillar_shear_ex(; highres = false)
  cyc = transit -> begin
    loc = almost_flat_caterpillar(transit)
    angle -> Caterpillar.cocycle(angle, loc)
  end
  perturbation = Matrix[
    traceless(1, 0, 0),
    traceless(0, 1, 0),
    traceless(0, 0, 1),
    traceless(-1, 0, 0)
  ]
  shear_ex(cyc, Regular.generators(4, 4), perturbation, highres = highres)
end

# === limit curves

function elem_generators(fp, a)
  last_fp = sum([fp[k] for k in 1:2:length(fp)]) - sum([fp[k] for k in 2:2:length(fp)])
  if iseven(length(fp))
    last_fp = -last_fp
  end
  [[a b; 0 1/a] for b in (1/a - a)*[fp; last_fp]]
end

function all_products(length, transit)
  if (length > 1)
    prev = all_products(length-1, transit)
    return vcat([[t*g for g in prev] for t in transit]...)
  else
    return transit
  end
end

# a dot at the most contracted eigenline of m
function eigen_dot(m, radius)
  m_eigen = eigen(m)
  x, y = reim(planeproj(m_eigen.vectors[:,argmin(map(abs, m_eigen.values))]))
  box = context(x - radius, y - radius, 2radius, 2radius)
  compose(box, Compose.circle())
end

# a dot at the most contracted line of m
function stable_dot(m, radius, color = nothing)
  x, y = reim(planeproj(stable(m)))
  box = context(x - radius, y - radius, 2radius, 2radius)
  if color == nothing
    compose(box, Compose.circle())
  else
    compose(box, Compose.circle(), fill(color))
  end
end

function elem_ex()
  # set parameters
  transit = elem_generators([cis(2*pi*k/6) for k in 0:4], sqrt(3))
  
  # draw limit set
  limit_set = [eigen_dot(m, 0.002) for m in all_products(5, transit)]
  
  # build and evolve cocycle
  in_lengths = [sin((@interval(k) - 1//5) * @interval(pi)/length(transit)) for k in 1:length(transit)]
  orig = Cocycle(cumsum(in_lengths), transit, collect(length(transit):-1:1))
  iter = power_twostep(orig, 8)
  
  # draw limit curve
  limit_curve = [
    stable_dot(bl.f_transit, 0.004, "tomato")
    for bl in iter.blocks_by_in
  ]
  
  # compose picture
  picture = compose(context(units = UnitBox(-1, -1, 2, 2)),
    (context(), limit_curve...),
    (context(), limit_set...),
    (context(), rectangle(), fill("white"))
  )
  draw(SVG("elem-limit.svg", 7cm, 7cm), picture)
end

# === abelianization map for local systems on a square

function square_ab_data(level_curves, cyc_family, window = false, neg = false)
  level_col = []
  x1_col = []
  x2_col = []
  
  flip = neg ? -1 : 1
  
  if window == false
    for i in 1:length(level_curves)
      for params in level_curves[i]
        exp_x = shears(cyc_family(params))
        push!(level_col, params[1])
        push!(x1_col, log(flip * real(exp_x[1])))
        push!(x2_col, log(flip * real(exp_x[2])))
      end
    end
  else
    for level in level_curves
      # march downward
      step = 0.005
      t = 1
      x = [0, 0]
      while true
        exp_x = shears(cyc_family((level, exp(t))))
        x = [log(flip * real(u)) for u in exp_x]
        if x[1] < 1.5*window[1] || 1.5*window[2] < x[1]
          step /= 2
          t += step
        else
          push!(level_col, level)
          push!(x1_col, x[1])
          push!(x2_col, x[2])
          if x[1] < 1.05*window[1] || 1.05*window[2] < x[1]
            break
          else
            t -= step
          end
        end
      end
      
      # march upward
      step = 0.005
      t = 1 + step
      x = [0, 0]
      while true
        exp_x = shears(cyc_family((level, exp(t))))
        x = [log(flip * real(u)) for u in exp_x]
        if x[1] < 1.5*window[1] || 1.5*window[2] < x[1] || x[2] < 1.5*window[3]
          step /= 2
          t -= step
        else
          pushfirst!(level_col, level)
          pushfirst!(x1_col, x[1])
          pushfirst!(x2_col, x[2])
          if x[1] < 1.05*window[1] || 1.05*window[2] < x[1] || x[2] < 1.05*window[3]
            break
          else
            t += step
          end
        end
      end
    end
  end
  
  DataFrame(
    Any[level_col, x1_col, x2_col],
    map(Symbol, ["level", "x1", "x2"])
  )
end

function square_ab_plot(data, angle, window, levelname)
  fade = t -> RGB(0/255, (30*(1-t) + 216*t)/255, (140*(1-t) + 180*t)/255)
  levels_layer = layer(data, x = :x1, y = :x2, color = :level, Geom.path)
  bdry_layer = layer(
    x = [1.05*window[1], 1.05*window[2]],
    y = [-t*mid(tan(angle)) for t in [1.05*window[1], 1.05*window[2]]],
    Geom.path,
    Theme(default_color = colorant"black")
  )
  plot(
    levels_layer, bdry_layer,
    Scale.color_asinh(colormap = fade),
    Coord.cartesian(
      xmin = window[1], xmax = window[2],
      ymin = window[3], ymax = window[4],
      fixed = true
    ),
    Guide.xlabel("log (<i>A</i><sub>ab</sub><sup>+</sup>)<sub>1</sub>"),
    Guide.ylabel("log (<i>A</i><sub>ab</sub><sup>+</sup>)<sub>2</sub>"),
    Guide.colorkey(title = levelname)
  )
end

function punkd_torus_plot(; angle = @interval(3//7), svg = false)
  level_curves = Any[
    [(exp(s), t) for t in range(-6.2, stop = 6.2, length = 64)]
    for s in range(-2, stop = 2, length = 30)
  ]
  fam = p -> begin
    loc = PunkdTorusLocSys(p[1], p[2])
    Square.cocycle(angle, loc)
  end
  window = (-3, 3, -3, 3)
  data = square_ab_data(level_curves, fam)
  plot = square_ab_plot(data, angle, window, "length")
  plot |> (svg ? SVG : PDF)("punkd-torus-ab-map.$(svg ? "svg" : "pdf")", 20cm, 16cm)
end

function quasicrystal_plot(; angle = @interval(3//7), svg = false)
  level_curves = vcat(range(0.4, stop = 4.4, length = 21), range(-0.4, stop = -1.4, length = 6))
  fam = p -> begin
    loc = QuasicrystalLocSys(p[1], p[2])
    Square.cocycle(angle, loc)
  end
  window = (-2, 2, -2, 2)
  data = square_ab_data(level_curves, fam, window, true)
  plot = square_ab_plot(data, angle, (-2, 1, -2, 2), "binding energy")
  plot |> (svg ? SVG : PDF)("quasicrystal-ab-map.$(svg ? "svg" : "pdf")", 20cm, 16cm)
end

# === geodesic lamination movie

# the first two terms of a sawtooth wave, modified to zero out the jerk at the
# the inflection point and rescaled into the box [0,1] × [0,1]
function easing(t)
  θ = π*(t - 1/2)
  x = sin(θ) - sin(3θ)/26
  (1 + (x / (1+1/26))) / 2
end

function movie(; shape = CaterpillarLocSys, ascent = 2, eps = 1e-3, theme = tacos, testframe = true, center = nothing, svg = false, verbose = true)
  if shape == CaterpillarLocSys
    # enumerate symmetry group elements
    transit = Regular.generators(4, 4)
    dbl_transit = [transit; [inv(t) for t in transit]]
    
    # set up local system
    loc = almost_flat_caterpillar(transit);
    
    # set up crawler
    crawler = TileCrawler(4, 4, ascent)
    findhome!(crawler, dbl_transit)
    
    # set angle range
    test_angle = @interval(3π/4 + 1//11)
    start = @interval(358//114)
    fin = @interval(180//114)
  elseif shape == PunkdTorusLocSys
    # set up local system
    loc = PunkdTorusLocSys(2, 1/3)
    
    # set up crawler
    crawler = FreeCrawler(2, 5)
    findhome!(crawler, [loc.n_transit, loc.e_transit, loc.s_transit, loc.w_transit])
    
    # set angle range
    test_angle = @interval(3//7)
    start = @interval(179//114)
    fin = @interval(1//114)
  end
  
  # render frames
  if testframe
    println("Test frame")
    frame = @time(Lamination.render(test_angle, loc, crawler, eps, theme, center = center, verbose = verbose))
    if svg
      draw(SVG("triangle_test.svg", 7cm, 7cm), frame)
    else
      draw(PDF("triangle_test.pdf", 7cm, 7cm), frame)
    end
  else
    n = 25
    for t in 0:n
      println("Frame $t")
      u = easing(@interval(t//n))
      frame = @time(Lamination.render(@interval((1-u)*start + u*fin), loc, crawler, eps, theme, center = center, verbose = verbose))
      draw(PNG(@sprintf("triangle_mov/frame%02i.png", t), 500px, 500px), frame)
    end
  end
end

# === latitude geodesic on a punctured torus

# draw some lifts of the geodesic around the latitude of a punctured torus. note
# that you can't hit all the lifts just by increasing the ascent of the crawler.
function latitude_geodesic(; eps = 1e-3, theme = shell, svg = false)
  # set up local system
  loc = PunkdTorusLocSys(2, 1/3)
  
  # set up crawler
  crawler = FreeCrawler(4, 3)
  branches = [loc.e_transit^k * loc.n_transit * loc.w_transit^k for k in -1:2]
  findhome!(crawler, vcat(branches, map(inv, branches)))
  
  # start a list of layers
  layers = []
  
  # draw geodesic lifts
  stable_lines = [stable(loc.e_transit), stable(loc.w_transit)]
  base = IdealPolygon(1, [planeproj(v) for v in stable_lines])
  lifts = mapcollect(orbiter(base, eps, polygon_penciller(theme), diam = [1, 2]), crawler)
  lift_gp = compose(context(), lifts..., stroke(theme.leafcolor), linewidth(0.1mm), fill(nothing))
  push!(layers, lift_gp)
  
  # draw background
  if theme.diskcolor != nothing
    disk = compose(context(), Compose.circle(), fill(theme.diskcolor), stroke(nothing))
    push!(layers, disk)
  end
  
  # render
  picture = compose(context(), layers...)
  if svg
    picture |> SVG("latitude_geodesic.svg", 7cm, 7cm)
  else
    picture |> PDF("latitude_geodesic.pdf", 7cm, 7cm)
  end
end

# === ideal triangulation of a punctured torus ===

const strawberry_sunrise = [
  RGB(255/255, 118/255, 188/255),
  RGB(255/255, 200/255, 122/255),
  RGB(245/255, 234/255, 215/255)
]

const chocolate = [
  RGB(53/255, 16/255, 3/255),
  RGB(118/255, 82/255, 62/255),
  RGB(161/255, 127/255, 103/255)
]

function torus_punks(h)
  p = möbius_map([[im, -1] [-1, im]], exp(h/2))
  return [p, -conj(p), -p, conj(p)]
end

# Yair Minsky, "The classification of punctured-torus groups"
# Equation 2.1
function torus_sym(h)
  l_frame = [[1, 1] [-1, 1]]
  h_frame = [[1im, 1] [-1im, 1]]
  l = 2atanh(sech(h/2))
  return Array[
    l_frame * Diagonal([exp(-l/2), exp(l/2)]) * inv(l_frame),
    h_frame * Diagonal([exp(-h/2), exp(h/2)]) * inv(h_frame)
  ]
end

# given a number n, return some functions that each take a möbius transformation
# m and apply it to part ideal triangulation of a punctured torus, Dehn-twisted
# n times in the `down` direction
function flip_orbiters(n, h, down)
  up = inv(down)
  
  # write down the initial ends of the edges
  p = torus_punks(h)
  x_start, x_end = p[4], p[3]
  y_start, y_end = p[4], p[2]
  z_start, z_end = p[4], p[1]
  
  # flip the edges
  is_x_turn = true
  for c in 1:n
    if is_x_turn
      x_start = möbius_map(up, x_start)
      x_end = möbius_map(down, x_end)
    else
      y_start = möbius_map(up, y_start)
      y_end = möbius_map(down, y_end)
    end
    is_x_turn = !is_x_turn
  end
  
  # record the vertices of the central lifts of the two triangles
  if iseven(n)
    q = [-x_end, y_end, x_end, y_start]
  else
    q = [x_start, y_start, -x_start, y_end]
  end
  
  # return
  [
    # x edge
    m -> begin
      shift_x_start = möbius_map(m, x_start)
      shift_x_end = möbius_map(m, x_end)
      compose(context(), ideal_edges(shift_x_start, shift_x_end))
    end,
    # y edge
    m -> begin
      shift_y_start = möbius_map(m, y_start)
      shift_y_end = möbius_map(m, y_end)
      compose(context(), ideal_edges(shift_y_start, shift_y_end))
    end,
    # z edge
    m -> begin
      shift_z_start = möbius_map(m, z_start)
      shift_z_end = möbius_map(m, z_end)
      compose(context(), ideal_edges(shift_z_start, shift_z_end))
    end,
    m -> begin
      shift_q = [möbius_map(m, u) for u in q]
      compose(
        context(),
        # dark triangle foliation
        (
          context(),
          horotriangle(shift_q[1], shift_q[4], shift_q[2], 69, 1/21, 4e-3),
          stroke(chocolate[2])
        ),
        # light triangle foliation
        (
          context(),
          horotriangle(shift_q[3], shift_q[2], shift_q[4], 69, 1/21, 4e-3),
          stroke(chocolate[3])
        )
      )
    end
  ]
end

function render_flip(crawler::CayleyCrawler, orbiters, name, foliate = false, svg = false)
  # draw background and boundary
  disk = compose(context(), Compose.circle(), fill(chocolate[1]), stroke(nothing))
  bdry = compose(context(), Compose.circle(), stroke("white"), linewidth(0.25mm), fill(nothing))
  
  # draw lamination
  lam_gps = []
  for i in 1:3
    lam_edges = mapcollect(orbiters[i], crawler)
    push!(lam_gps, compose(
      context(),
      lam_edges...,
      stroke(strawberry_sunrise[i]),
      linewidth(0.3mm),
      fill(nothing)
    ))
  end
  
  # draw foliation
  if foliate
    fol = mapcollect(orbiters[4], crawler)
    fol_gp = compose(context(), fol..., linewidth(0.25mm))
    picture = compose(context(), bdry, (lam_gps..., fol_gp), disk)
  else
    picture = compose(context(), bdry, (lam_gps...,), disk)
  end
  
  # render
  if svg
    draw(SVG(name, 9cm, 9cm), picture)
  else
    draw(PDF(name, 9cm, 9cm), picture)
  end
end

function various_shears()
  crawler = FreeCrawler(2, 5)
  h = range(2, stop = 4, length = 5)
  for i in 1:length(h)
    # set up holonomy group crawler
    left, down = torus_sym(h[i])
    dbl_transit = [left, down, inv(left), inv(down)]
    findhome!(crawler, dbl_transit)
    
    # render
    if i == 1
      render_flip(
      crawler,
      flip_orbiters(0, h[i], down),
      "shears/disk1.pdf",
      false
    )
    render_flip(
      crawler,
      flip_orbiters(0, h[i], down),
      @sprintf("shears/shear%i.pdf", i),
      true
    )
    end
  end
end

function animate_flips()
  h = 2.0
  
  # set up holonomy group crawler
  crawler = FreeCrawler(2, 5)
  left, down = torus_sym(h)
  dbl_transit = [left, down, inv(left), inv(down)]
  findhome!(crawler, dbl_transit)
  
  # render flips
  for frame in 0:6
    render_flip(
      crawler,
      flip_orbiters(frame, h, down),
      @sprintf("flips/flip%i.pdf", frame+1)
    )
  end
end

end

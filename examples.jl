include("interval_exchange.jl")
include("caterpillar.jl")
include("rectangle.jl")
include("regular.jl")
include("cayley_crawler.jl")

module Examples

using
  Gadfly,
  DataFrames,
  Compose,
  Colors,
  ValidatedNumerics,
  IntervalExchange,
  Rectangle,
  Caterpillar,
  PoincaréDisk,
  Crawl

import Regular

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

# print a "twisted caterpillar" cocycle and its abelianization
function abelianization_ex()
  # set up cocycle
  orig = twisted_caterpillar(@interval(3π/4 + 1//11), Regular.generators(4, 4))
  
  # evolve cocycle
  iter = orig
  for _ in 1:4
    iter = twostep(iter)
  end
  
  # abelianize cocycle
  ab = abelianize(orig, iter)
  
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
  iter = orig
  for _ in 1:4
    iter = twostep(iter)
  end
  
  # abelianize
  ab = abelianize(orig, iter)
  [bl.f_transit[1,1] for bl in ab.blocks_by_in]
end

# linspace doesn't work with Interval objects, so here's a slapdash replacement
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
        append!(angle_col, collect(repeated(mid(angles[i]), length(x))))
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
    default_color = tacos[1],
    default_point_size = 0.5mm,
    highlight_width = 0mm,
  )
  imag_theme = Theme(
    default_color = tacos[4],
    default_point_size = 0.5mm,
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
    expm(t*a)*g*expm(-t*a)
  end

function shear_ex(cyc, transit, perturbation; highres = false)
  # no perturbation
  println("=== no perturbation\n")
  data_no = shear_data(cyc(transit), highres ? 300 : 18)
  p_no = shear_plot(data_no)
  page_height = length(transit)*3cm + 1cm
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
    loc = RectangleLocSys{AbstractInterval}(
      @interval(1), @interval(1),
      inv(transit[2]), transit[1]
    )
    angle -> cocycle(angle, loc)
  end
  perturbation = Matrix[
    traceless(10, 0, 0),
    traceless(0, 10, 0)
  ]
  shear_ex(cyc, Regular.generators(2, k), perturbation, highres = highres)
end

function caterpillar_shear_ex(; highres = false)
  cyc = transit -> (angle -> twisted_caterpillar(angle, transit))
  perturbation = Matrix[
    traceless(1, 0, 0),
    traceless(0, 1, 0),
    traceless(0, 0, 1),
    traceless(-1, 0, 0)
  ]
  shear_ex(cyc, Regular.generators(4, 4), perturbation, highres = highres)
end

# === geodesic lamination movie

# color scheme
const tacos = [
  RGB(255/255, 1/255, 73/255),
  RGB(255/255, 121/255, 1/255),
  RGB(255/255, 210/255, 0/255),
  RGB(0/255, 200/255, 146/255)
]

# the first two terms of a sawtooth wave, modified to zero out the jerk at the
# the inflection point and rescaled into the box [0,1] × [0,1]
function easing(t)
  θ = π*(t - 1/2)
  x = sin(θ) - sin(3θ)/26
  (1 + (x / (1+1/26))) / 2
end

# given a jump j, return the function that takes a möbius transformation m and
# applies it to the triangle associated with j
orbiter(j::Jump) =
  m -> begin
    left = möbius_map(m, planeproj(j.left_stable))
    right = möbius_map(m, planeproj(j.right_stable))
    pivot = möbius_map(m, planeproj(j.pivot_stable))
    compose(
      context(),
      (context(),
        ideal_path(left, right, pivot),
        fill(tacos[j.sing])
      )
    )
  end

# draw the complementary triangles of the geodesic lamination specified by the
# "twisted caterpillar" cocycle with the given angle and transition maps. if a
# frame number is specified, render an appropriately named bitmap to be used as
# a frame of a movie. otherwise, render a PDF test frame
function render{R <: AbstractInterval}(
  angle::R,
  transit,
  crawler::CayleyCrawler,
  orbiter;
  frame = nothing
)
  # set up cocycle
  orig = twisted_caterpillar(angle, transit)
  
  # evolve cocycle
  iter = orig
  for i in 1:4
    print("  Step $i\n  ")
    iter = @time(twostep(iter))
    println("    $(length(iter.blocks_by_in)) blocks")
  end
  
  # find the widest triangle for each singularity
  b_jumps = scancollect(iter, Jump, b_fn = BJump)
  widest = Jump[]
  for sing in 1:4
    push!(widest, maximum(filter(j -> j.sing == sing, b_jumps)))
  end
  
  # draw background
  bg = compose(context(),
    (context(), rectangle(), fill("white"), stroke(nothing))
  )
  
  # draw triangle lifts
  triangles = vcat([mapcollect(orbiter(j), crawler) for j in widest]...)
  lam_cmp = compose(context(), triangles...)
  
  # render
  picture = compose(lam_cmp, bg)
  if frame == nothing
    draw(PDF("triangle_test.pdf", 7cm, 7cm), picture)
  else
    draw(PNG(@sprintf("triangle_mov/frame%02i.png", frame), 500px, 500px), picture)
  end
end

function movie(; testframe = true)
  # enumerate symmetry group elements
  transit = Regular.generators(4, 4)
  dbl_transit = [transit; [inv(t) for t in transit]]
  crawler = TileCrawler(4, 4, 2)
  findhome!(crawler, dbl_transit)
  
  if testframe
    println("Test frame")
    @time(render(@interval(3π/4 + 1//11), transit, crawler, orbiter))
  else
    start = @interval(358//114)
    fin = @interval(180//114)
    n = 25
    for t in 0:n
      println("Frame $t")
      u = easing(@interval(t//n))
      @time(render(@interval((1-u)*start + u*fin), transit, crawler, orbiter, frame = t))
    end
  end
end

# === ideal triangulation of a punctured torus ===

const strawberry_sunrise = [
  RGB(70/255, 23/255, 7/255),
  RGB(255/255, 118/255, 188/255),
  RGB(255/255, 200/255, 122/255),
  RGB(245/255, 234/255, 215/255)
]

# given a number n, return the function that takes a möbius transformation m and
# applies it to the an ideal triangulation of a punctured torus, Dehn-twisted n
# times in the `down` direction
function lam_orbiter(n, down)
  vertex = [cis((2j+1)*π/4) for j in 0:3]
  for c in 1:n
    vertex[2] = möbius_map(down, vertex[2])
    vertex[3] = möbius_map(down, vertex[3])
  end
  for c in 1:n
    vertex[1] = möbius_map(inv(down), vertex[1])
    vertex[4] = möbius_map(inv(down), vertex[4])
  end
  m -> begin
    w = [möbius_map(m, v) for v in vertex]
    compose(
      context(),
      (context(), ideal_edges(w[1], w[2]), stroke(strawberry_sunrise[2])),
      (context(), ideal_edges(w[1], w[3]), stroke(strawberry_sunrise[3])),
      (context(), ideal_edges(w[1], w[4]), stroke(strawberry_sunrise[4]))
    )
  end
end

function triangulate()
  down, right = Regular.generators(2, nothing)
  dbl_transit = [down, right, inv(down), inv(right)]
  crawler = FreeCrawler(2, 6)
  findhome!(crawler, dbl_transit)
  
  # draw background and boundary
  bg = compose(context(), rectangle(), fill("white"), stroke(nothing))
  disk = compose(context(), circle(), fill(strawberry_sunrise[1]), stroke(nothing))
  bdry = compose(context(), circle(), stroke("white"), linewidth(0.4mm), fill(nothing))
  
  # draw lamination
  lam_edges = mapcollect(lam_orbiter(1, down), crawler)
  lam_cmp = compose(context(), lam_edges..., linewidth(0.2mm), fill(nothing))
  
  # render
  lam_picture = compose(context(), (context(1/9, 1/9, 7/9, 7/9), bdry, lam_cmp, disk), bg)
  ##fol_picture = compose(context(), (context(1/9, 1/9, 7/9, 7/9), bdry, lam_cmp, fol_cmp), bg)
  draw(SVG("laminated.svg", 9cm, 9cm), lam_picture)
  ##draw(SVG("foliated.svg", 9cm, 9cm), fol_picture)
end

end

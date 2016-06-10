include("caterpillar.jl")
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
  orig = twisted_caterpillar(@interval(3π/4 + 1//11), Regular.generators(2))
  
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

# compute shear parameters of a "twisted caterpillar" cocycle, with the given
# transition maps, over a range of angles, with the given resolution. the range
# of angles goes roughly from π/2 to π. we use rational approximations for π/2
# and π that are only accurate to about one part in 114, making it unlikely that
# we'll get within machine precision of a saddle connection at any reasonable
# resolution.
function shear_data(transit, res::Integer)
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
        
        x = shears(twisted_caterpillar(angles[i], transit))
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

function shear_plot_ex(; highres = false)
  perturbation = Matrix[
    traceless(1, 0, 0),
    traceless(0, 1, 0),
    traceless(0, 0, 1),
    traceless(-1, 0, 0)
  ]
  
  # no perturbation
  println("=== no perturbation\n")
  data_no = shear_data(Regular.generators(2), highres ? 300 : 18)
  p_no = shear_plot(data_no)
  draw(PDF("no-perturbation.pdf", 30cm, 40cm), p_no)
  println()
  
  # small perturbation
  println("=== small perturbation\n")
  transit_sm = map(expconj(0.01), zip(Regular.generators(2), perturbation))
  data_sm = shear_data(transit_sm, highres ? 300 : 18)
  p_sm = shear_plot(data_sm)
  draw(PDF("small-perturbation.pdf", 30cm, 40cm), p_sm)
  println()
  
  # large perturbation
  println("=== large perturbation\n")
  transit_lg = map(expconj(0.1), zip(Regular.generators(2), perturbation))
  data_lg = shear_data(transit_lg, highres ? 300 : 18)
  p_lg = shear_plot(data_lg)
  draw(PDF("large-perturbation.pdf", 30cm, 40cm), p_lg)
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
  transit = Regular.generators(2)
  dbl_transit = [transit; [inv(t) for t in transit]]
  crawler = CayleyCrawler(4, 4, 2)
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

end

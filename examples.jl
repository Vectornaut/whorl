include("regular.jl")
include("cayley_crawler.jl")
include("caterpillar.jl")

module Examples

using Gadfly, Compose, Colors, ValidatedNumerics, IntervalExchange, Caterpillar, PoincaréDisk, Regular, Crawl

# === abelianized holonomy plots

function ab_hols(orig)
  # evolve cocycle
  iter = orig
  for i in 1:4
    iter = twostep(iter)
  end
  
  # abelianize
  ab = abelianize(orig, iter)
  [real(bl.f_transit[1,1]) for bl in ab.blocks_by_in]
end

make_sl(a::Number, b::Number, c::Number) = [a b; c (1 + b*c)/a]

function remix(hols)
  new_hols = [
    hols[6],
    hols[8]*hols[6],
    hols[6]*hols[8]*hols[13]*hols[2]*hols[12]*hols[3]*hols[11]*hols[5]*hols[9],
    hols[8]^3*hols[6]^3*hols[1]^2*hols[5]*hols[7]*hols[12]*hols[13]
  ]
end

# linspace doesn't work with Interval objects, so here's a slapdash replacement
function grid(start, fin, res::Integer)
  map(u -> (1-u)*start + u*fin, [@interval(t//res) for t in 0:res])
end

function regular_table(res::Integer)
  transit = generators(2)
  domain = grid(@interval(358//114), @interval(180//114), res)
  hol_list = []
  for angle in domain
    print("=")
    push!(hol_list, ab_hols(twisted_caterpillar(angle, transit)))
  end
  println()
  (domain, hcat(hol_list...))
end

function poster_table(res::Integer)
  transit = Matrix[
    make_sl(2.1, 0.1, 0.2),
    make_sl(1.9, 0.1, -0.1),
    make_sl(2.2, -0.1, 0),
    make_sl(1.8, 0.1, -0.2),
    make_sl(2.3, 0, 0.2),
    make_sl(0.8, 0.7, -0.6),
    make_sl(0.9, -0.5, 0.7)
  ]
  domain = grid(@interval(358//114), @interval(180//114), res)
  hol_list = []
  for angle in domain
    print("=")
    push!(hol_list, remix(ab_hols(caterpillar_cocycle(angle, transit...))))
  end
  println()
  (domain, hcat(hol_list...))
end

function poster_plot(domain, hol_table, hsize, vsize, name; special_ranges = false)
  yranges = [0 1.3; 0 0.5; -0.05 0.01; -0.01 0.02]
  ranges = [
    Coord.Cartesian(
      xmin = Float64(π/2), xmax = Float64(π),
      ymin = yranges[i,1], ymax = yranges[i,2]
    ) for i in 1:size(yranges, 1)
  ]
  graphs = [
    plot(
      x = map(mid, domain),
      y = hol_table[i,:],
      special_ranges ? ranges[i] : Coord.Cartesian(),
      Geom.point,
      Guide.xlabel(nothing),
      Guide.ylabel("Coordinate $i"),
      Theme(
        default_color = parse(Colorant, "black"),
        default_point_size = 0.75mm,
        highlight_width = 0mm,
        major_label_font = "Anaheim"
      ),
      Guide.xticks(
        ticks = [
          map(
            t -> atan(t) + π/2,
            [0, 1/4, 1/3, 1/2, 2/3, 3/4, 1, 4/3, 3/2, 2/1, 3/1, 4/1]
          )
          π
        ],
        label = false
      ),
      special_ranges ? Guide.yticks(ticks = yranges[i,:]) : Guide.yticks()
    ) for i in 1:size(hol_table, 1)
  ]
  draw(PDF(name, hsize, vsize), vstack(graphs...))
end

function regular_example()
  domain, hol_table = regular_table(300)
  poster_plot(domain, hol_table, 44cm, 66cm, "regular-example.pdf")
end

function poster_example()
  domain, hol_table = poster_table(300)
  poster_plot(domain, hol_table, 44cm, 22cm, "poster-example.pdf", special_ranges = true)
end

# === geodesic lamination movie

tacos = [
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

function render{R <: AbstractInterval}(angle::R, transit, crawler::CayleyCrawler, orbiter; frame = nothing)
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
  transit = generators(2)
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

using Colors

include("poincare_disk.jl")
include("color_scheme.jl")

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
    l_frame * diagm([exp(-l/2), exp(l/2)]) * inv(l_frame),
    h_frame * diagm([exp(-h/2), exp(h/2)]) * inv(h_frame)
  ]
end

function guides(h_old, h_new)
  p = torus_punks(h_old)
  q = torus_punks(h_new)
  
  # write down the möbius transformation
  #  ∞ --> p[2]
  #  0 --> p[3]
  #  1 --> p[1]
  m2 = [[p[2]*(p[1] - p[3]), p[1] - p[3]] [p[3]*(p[2] - p[1]), p[2] - p[1]]]
  
  # write down the möbius transformation
  #  ∞ --> q[2]
  #  0 --> q[3]
  #  1 --> q[1]
  n2 = [[q[2]*(q[1] - q[3]), q[1] - q[3]] [q[3]*(q[2] - q[1]), q[2] - q[1]]]
  
  # write down the möbius transformation
  #  ∞ --> p[4]
  #  0 --> p[3]
  #  1 --> p[1]
  m4 = [[p[4]*(p[1] - p[3]), p[1] - p[3]] [p[3]*(p[4] - p[1]), p[4] - p[1]]]
  
  # write down the möbius transformation
  #  ∞ --> q[4]
  #  0 --> q[3]
  #  1 --> q[1]
  n4 = [[q[4]*(q[1] - q[3]), q[1] - q[3]] [q[3]*(q[4] - q[1]), q[4] - q[1]]]
  
  # apply the mobius transformations to get the desired guide lines
  return compose(
    context(),
    geodesic(q[2], möbius_map(n2*inv(m2), p[4])),
    geodesic(q[2], möbius_map(n2*inv(m2), cis(π/7)*p[4])),
    geodesic(q[2], möbius_map(n2*inv(m2), cis(-π/7)*p[4])),
    geodesic(q[4], möbius_map(n4*inv(m4), p[2])),
    geodesic(q[4], möbius_map(n4*inv(m4), cis(π/7)*p[2])),
    geodesic(q[4], möbius_map(n4*inv(m4), cis(-π/7)*p[2]))
  )
end

function torus_edges(h, depth)
  p = torus_punks(h)
  sym = torus_sym(h)
  
  return compose(
    context(),
    geodesic_orbit(p[1], p[4], sym, depth),
    geodesic_orbit(p[1], p[3], sym, depth),
    geodesic_orbit(p[1], p[2], sym, depth)
  )
end

function torus_chart(h)
  p = torus_punks(h)
  sym = torus_sym(h)
  
  return compose(
    context(),
    geodesic(p[3], p[4]),
    geodesic(p[1], p[4]),
    geodesic(p[1], p[3]),
    geodesic(p[1], p[2]),
    geodesic(p[3], p[2])
  )
end

function torus_foliage(h, depth)
  p = torus_punks(h)
  sym = torus_sym(h)
  
  return compose(
    context(),
    horotriangle(p[4], p[3], p[1], 56, 4/49, 4e-3),
    horotriangle(p[2], p[1], p[3], 56, 4/49, 4e-3)
  )
end

function triangle_foliage()
  p = [im*cis(2π/3*k) for k in 1:3]
  return compose(
    context(),
    compose(
      context(),
      geodesic(p[1], p[2]),
      geodesic(p[2], p[3]),
      geodesic(p[3], p[1]),
      stroke(vert_ink)
    ),
    compose(
      context(),
      horotriangle(p[1], p[2], p[3], 56, 4/49, 4e-3),
      stroke(hor_ink)
    ),
    linewidth(0.75pt)
  )
end

function torus_pics()
  h_old = 3
  h_new = 2
  h_fol = 2.6
  
  # draw poincaré disk
  disk = compose(context(), circle(), fill(disk_ink), linewidth(0.4mm))
  
  # draw lamination and foliation
  lam_orig = compose(context(), torus_edges(h_old, 5), stroke(vert_ink), linewidth(0.75pt))
  lam_shift = compose(context(), torus_edges(h_new, 5), stroke(vert_ink), linewidth(0.75pt))
  lam_fol = compose(context(), torus_edges(h_fol, 5), stroke(vert_ink), linewidth(0.75pt))
  fol = compose(context(), torus_foliage(h_fol, 5), stroke(hor_ink), linewidth(0.75pt))
  
  # draw torus chart
  chart = compose(context(), torus_chart(h_new), stroke(vert_ink), linewidth(0.75pt))
  
  # draw guides
  gl = compose(context(), guides(h_old, h_old), stroke("red"), linewidth(0.25pt))
  gl_shift = compose(context(), guides(h_old, h_new), stroke("red"), linewidth(0.25pt))
  
  # print outputs
  lam_file = SVG("punk-lam-orig_raw.svg", 5.2cm, 5.2cm)
  lam_shift_file = SVG("punk-lam-shifted_raw.svg", 5.2cm, 5.2cm)
  fol_file = SVG("punk-fol_raw.svg", 5.2cm, 5.2cm)
  chart_file = SVG("chart-shifted_raw.svg", 5.2cm, 5.2cm)
  tri_fol_file = SVG("tri-fol_raw.svg", 5.2cm, 5.2cm)
  draw(lam_file, compose(context(), gl, lam_orig, disk))
  draw(lam_shift_file, compose(context(), gl_shift, lam_shift, disk))
  draw(fol_file, compose(context(), lam_fol, fol, disk))
  draw(tri_fol_file, compose(context(), triangle_foliage(), disk))
  draw(chart_file, compose(context(), gl_shift, chart, disk))
end

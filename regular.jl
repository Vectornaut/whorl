module Regular

export generators

using LinearAlgebra

##[type cleanup] use a concrete type like Int instead of Integer for counting to
## improve performance

# the horizontal translation of the Poincaré disk that preserves the tesselation
# with Schläfli symbol {2j, 2k}. for k, a value of nothing is treated as
# infinity.
function slide(j::Integer, k::Union{Integer, Nothing}, uhp)
  # the hyperbolic cosine of the inscribed radius of a face, from the formula in
  # section 5 of Coxeter's article on "the trigonometry of hyperbolic
  # tesselations"
  ##cosh_r = (k == nothing ? 1 : cos(pi/(2k))) / sin(pi/(2j))
  cosh_r = (k == nothing ? 1 : cos(BigFloat(pi)/(2k))) / sin(BigFloat(pi)/(2j)) ##[higher precision]
  
  # the exponential of the inscribed radius of a face
  exp_r = cosh_r + sqrt(cosh_r^2 - 1)
  
  # the inscribed radius of a face happens to be half the translation length of
  # the slide that sends the right edge of the central face to the left edge.
  # we start with a transformation that realizes the slide in any half-plane
  # model, and then conjugate by a transformation that sends the left half-plane
  # to the unit disk if needed.
  frame = uhp ? I : [[1, 1] [-1, 1]]
  frame * Diagonal([1/exp_r, exp_r]) * inv(frame)
end

# a turn of s*π/j around the center of the Poincaré disk
function turn(j::Integer, s::Integer, uhp)
  if uhp
    [[cos(s*π/2j), sin(s*π/2j)] [-sin(s*π/2j), cos(s*π/2j)]]
  else
    Diagonal([cis(s*π/2j), cis(-s*π/2j)])
  end
end

# symmetry group generators for the tesselation with Schläfli symbol {2j, 2k}.
# for k, a value of nothing is treated as infinity.
generators(j::Integer, k::Union{Integer, Nothing}, uhp = false) =
  [turn(j, s, uhp) * slide(j, k, uhp) * turn(j, -s, uhp) for s in 1:j]

end

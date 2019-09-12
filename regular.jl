module Regular

export generators

using LinearAlgebra

##[type cleanup] use a concrete type like Int instead of Integer for counting to
## improve performance

# the horizontal translation of the Poincaré disk that preserves the tesselation
# with Schläfli symbol {2j, 2k}. for k, a value of nothing is treated as
# infinity.
function slide(j::Integer, k::Union{Integer, Nothing})
  # the hyperbolic cosine of the inscribed radius of a face, from the formula in
  # section 5 of Coxeter's article on "the trigonometry of hyperbolic
  # tesselations"
  cosh_r = (k == nothing ? 1 : cos(pi/(2k))) / sin(pi/(2j))
  
  # the exponential of the inscribed radius of a face
  exp_r = cosh_r + sqrt(cosh_r^2 - 1)
  
  # the inscribed radius of a face happens to be half the translation length of
  # the slide that sends the right edge of the central face to the left edge
  frame = [[1, 1] [-1, 1]]
  frame * diagm([1/exp_r, exp_r]) * inv(frame)
end

# a turn of s*π/j around the center of the Poincaré disk
turn(j::Integer, s::Integer) = [[cis(s*π/2j), 0] [0, cis(-s*π/2j)]]

# symmetry group generators for the tesselation with Schläfli symbol {2j, 2k}.
# for k, a value of nothing is treated as infinity.
generators(j::Integer, k::Union{Integer, Nothing}) = [turn(j, s) * slide(j, k) * turn(j, -s) for s in 1:j]

end

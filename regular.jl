module Regular

export generators

# the horizontal translation of the Poincaré disk that preserves the tesselation
function slide(genus::Integer)
  u = (cos(π/4genus) - sin(π/4genus)) / sqrt(cos(π/2genus)) # from laura's notebook
  λ = (1 - u)/(1 + u) # the translation length
  frame = [[1, 1] [-1, 1]]
  frame * diagm([λ, 1/λ]) * inv(frame)
end

turn(genus::Integer, s::Integer) = [[cis(s*π/2genus), 0] [0, 1]]

generators(genus::Integer) = [turn(genus, s) * slide(genus) * turn(genus, -s) for s in 1:2genus]

end

# the horizontal translation of the Poincaré disk that preserves the tesselation
function slide(genus::Integer)
  u = (cos(π/(4genus)) - sin(π/(4genus))) / sqrt(cos(π/(2genus))) # from laura's notebook
  λ = (1 - u)/(1 + u) # the translation length
  frame = [[1, 1] [-1, 1]]
  frame * diagm([λ, 1/λ]) * inv(frame)
end

turn(genus::Integer, k::Integer) = [[exp(1im*k*π/(2*genus)), 0] [0, 1]]

generators(genus::Integer) = [turn(genus, k) * slide(genus) * turn(genus, -k) for k in 1:2genus]

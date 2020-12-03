module Quasicrystal

using Main.Square

export QuasicrystalLocSys

struct QuasicrystalLocSys <: SquareLocSys
  n_transit
  e_transit
  s_transit
  w_transit
  
  # following Bjerklöv in "Dynamics of the Quasi-Periodic Schrödinger Cocycle",
  # we make the hopping matrix elements negative, as you'd get in a crystal with
  # s-type valence orbitals (see Nakano's notes, "Tight-Binding Model of
  # Electronic Structures"). the potential is zero at east-labeled sites and the
  # negative of `binding` at north-labeled sites.
  function QuasicrystalLocSys(binding, level; disk = false)
    transit = [
      [[-binding - level, 1] [-1, 0]],
      [[-level, 1] [-1, 0]],
      [[0, -1] [1, -binding - level]],
      [[0, -1] [1, -level]]
    ]
    if disk
      # the mobius transformation
      # -i --> -1
      #  1 -->  0
      #  i -->  1
      flatten = [[-im, 1] [im, 1]]
      
      transit = [inv(flatten) * g * flatten for g in transit]
    end
    new(transit[1], transit[2], transit[3], transit[4])
  end
end

end

module Quasicrystal

using Square

export QuasicrystalLocSys

type QuasicrystalLocSys <: SquareLocSys
  n_transit
  e_transit
  s_transit
  w_transit
  
  # following Bjerklöv in "Dynamics of the Quasi-Periodic Schrödinger Cocycle",
  # we make the hopping matrix elements negative, as you'd get in a crystal with
  # s-type valence orbitals (see Nakano's notes, "Tight-Binding Model of
  # Electronic Structures"). the potential is zero at east-labeled sites and the
  # negative of `binding` at north-labeled sites.
  QuasicrystalLocSys(binding, level) = new(
    [[-binding - level, 1] [-1, 0]],
    [[-level, 1] [-1, 0]],
    [[0, -1] [1, -binding - level]],
    [[0, -1] [1, -level]]
  )
end

end

module PunkdTorus

export PunkdTorusLocSys

using LinearAlgebra, Main.Square

struct PunkdTorusLocSys <: SquareLocSys
  # transitions
  n_transit
  e_transit
  s_transit
  w_transit
  
  # fundamental domain vertices, clockwise from northwest
  punks
  
  function PunkdTorusLocSys(length, twist)
    # the mobius transformation
    # -i --> -1
    #  1 -->  0
    #  i -->  1
    flatten = [[-im, 1] [im, 1]]
    
    # see Parker and Parkkonen's "Coordinates for Quasi-Fuchsian Punctured Torus
    # Space", Figure 1.2. our `arcleft` is the inverse of their S', and our
    # `expand` is their T.
    arcleft = [
       cosh(length/2)       -cosh(length/2) + 1;
      -cosh(length/2) - 1    cosh(length/2)
    ]
    expand = [
       cosh(twist/2)*coth(length/4)   -sinh(twist/2);
      -sinh(twist/2)                   cosh(twist/2)*tanh(length/4)
    ]
    
    # conjugate
    n_transit = inv(flatten) * arcleft * flatten
    e_transit = inv(flatten) * expand * flatten
    s_transit = inv(flatten) * inv(arcleft) * flatten
    w_transit = inv(flatten) * inv(expand) * flatten
    
    # return
    new(
      n_transit,
      e_transit,
      s_transit,
      w_transit,
      [eigvecs(p)[:,1] for p in [
        s_transit * e_transit * n_transit * w_transit,
        w_transit * s_transit * e_transit * n_transit,
        n_transit * w_transit * s_transit * e_transit,
        e_transit * n_transit * w_transit * s_transit
      ]]
    )
  end
end

end

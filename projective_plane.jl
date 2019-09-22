module ProjectivePlane

export
  veronese_emb,
  möbius_map,
  ##pts_to_pts,
  ##stable, # already in PoincaréDisk
  ##planeproj,
  geodesic,
  ideal_edges,
  ideal_path,
  horotriangle,
  geodesic_orbit

using Compose, LinearAlgebra, Main.PoincaréDisk
using GenericLinearAlgebra ##[higher precision]

# === veronese embedding

function veronese_emb(m)
  # find special linear SVD. since m is special linear, the product of its
  # eigenvalues has to be 1. hence, no matter which U and Vt the SVD routine
  # finds, we know their determinants will have the same sign.
  svd_m = svd(m, full = true)
  if det(svd_m.U) > 0
    U = svd_m.U
    Vt = svd_m.Vt
  else
    U = svd_m.U[[2,1],:]
    Vt = svd_m.Vt[:,[2,1]]
  end
  
  xyt_to_lry = [-1 0 1; 1 0 1; 0 1 0]
  pre_rot = hcat(vcat(Vt^2, [0 0]), [0, 0, 1])
  boost = Diagonal([svd_m.S[1]^2, svd_m.S[2]^2, 1])
  post_rot = hcat(vcat(U^2, [0 0]), [0, 0, 1])
  post_rot * inv(xyt_to_lry) * boost * xyt_to_lry * pre_rot
end

# === projective transformations

##[type cleanup] switch to fixed arrays!

# apply a projective transformation, given as an operator on R^3, to a point on
# R^2.
## misnamed, obviously!
PoincaréDisk.möbius_map(m::Union{AbstractMatrix, UniformScaling}, v::Vector{T}) where T <: Number =
  [m[i,1]*v[1] + m[i,2]*v[2] + m[i,3] for i in 1:2] / (m[3,1]*v[1] + m[3,2]*v[2] + m[3,3])

# === points and geodesics
## misnamed: lines aren't geodesics in this case!

# project a line in affine space to the projective plane
planeproj(v::Vector{T}) where T <: Number = v[1:2] / v[3]

# a version of reim for use in Compose paths(z::Number) = (real(z)*cx, imag(z)*cy)
# turn a coordinate vector into a measure tuple
coord_measures(v::Vector{T}) where T <: Number = (v[1]*cx, v[2]*cy)

# draw the line between two points on the boundary of the projective plane
geodesic(tail::Vector{T}, head::Vector{T}) where T <: Number =
  compose(
    context(units=UnitBox(-1, -1, 2, 2)),
    line(coord_measures(tail), coord_measures(head))
  )

# draw an ideal polygon as a sequence of lines, good for stroking
## somewhat geometry-independent. put with more general code?
function PoincaréDisk.ideal_edges(verts::Vector{Vector{T}}) where T <: Number
  n = length(verts)
  cyc = i -> mod(i, n) + 1
  compose([geodesic(verts[cyc(i)], verts[cyc(i+1)]) for i in 0:(n-1)]...)
end

# draw an ideal polygon as a path, good for filling
## somewhat geometry-independent. put with more general code?
function PoincaréDisk.ideal_path(verts::Vector{Vector{T}}) where T <: Number
  n = length(verts)
  cyc = i -> mod(i, n) + 1
  compose(
    context(units=UnitBox(-1, -1, 2, 2)),
    polygon(coord_measures.(verts))
  )
end

end

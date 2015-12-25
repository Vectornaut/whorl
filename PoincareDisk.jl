module PoincareDisk
export geodesic

using Compose

# draws the geodesic between two points on the unit circle
function geodesic(tail::Complex, head::Complex)
  arcangle = sqrt(head/tail)
  radius = abs(imag(arcangle)/real(arcangle))
  
  if real(arcangle) > 1e-6
    # the geodesic is visibly curved
    arc = path([
            :M, real(tail)*cx, imag(tail)*cy,
            :A, radius*cx, radius*cy, 0, false, imag(arcangle) < 0, real(head)*cx, imag(head)*cy
          ])
  else
    # the geodesic looks straight
    arc = line([(real(tail)*cx, imag(tail)*cy), (real(head)*cx, imag(head)*cy)])
  end
  
  compose(context(units=UnitBox(-1, -1, 2, 2)), fill(nothing), arc)
end

draw(
  SVG("geode_test.svg", 10cm, 10cm),
  compose(
    context(),
    (context(), circle(), stroke("black"), fill(nothing)),
    (context(), geodesic(exp((π/5)im), exp(π*im)), stroke("yellowgreen")),
    (context(), geodesic(exp(0im), exp(π*im)), stroke("hotpink")),
  )
)

end # module

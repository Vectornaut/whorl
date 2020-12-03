using JLD

## list of examples
const HOPS_SHEARS = 1
const HOPS_SHEARS_ZOOM = 2
const SPIRALING_HOPS = 3
const GOLDEN_HOPS = 4

const files = [
  "hops/hops_shears.jld",
  "hops/hops_shears_zoom.jld",
  "hops/spiraling_hops.jld",
  "hops/golden_hops.jld"
]

const oldfiles = [
  "hops/old_hops_shears.jld",
  "hops/old_hops_shears_zoom.jld",
  "hops/old_spiraling_hops.jld",
  "hops/old_golden_hops.jld"
]

# https://github.com/JuliaLang/julia/blob/1303dfb96aa4345dca5d225e7ae72348307aef25/base/range.jl
struct LinSpace{T}
  start::T
  stop::T
  len::T
  divisor::T
end

# https://github.com/JuliaLang/julia/pull/25896
JLD.readas(a::LinSpace) = range(a.start, stop = a.stop, length = Int(a.len))

# https://github.com/JuliaIO/JLD.jl/blob/master/doc/jld.md
# commit d1cbe4a
translate("Base.LinSpace{Core.Float64}", "LinSpace{Core.Float64}")

function updatefile(example)
  data = load(oldfiles[example])
  save(
    files[example],
    "e_list", data["e_list"],
    "l_list", data["l_list"],
    "x_grid", data["x_grid"]
  )
end

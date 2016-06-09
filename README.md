# whorl

Whorl is software for playing with cocycles over interval exchange transformations. You can use it to do fun computations in geometric topology, including:

* Drawing pictures of geodesic laminations.
* Calculating shear coordinates of compact hyperbolic surfaces.
* [Abelianizing](https://www.ma.utexas.edu/users/afenyes/writing.html) local systems on compact translation surfaces.

I think Whorl has become useful and usable enough to be worth releasing, but it's still very immature, and future revisions are likely to involve major changes to its organization. You're welcome to check it out and play with it, and you could even use it to make figures or do calculations, but pay attention to which revision you're using&#x2014;pulling updates down from here could break your code without warning.

## Setup

1. Install [Julia](http://julialang.org/), the open-source technical computing language Whorl is written in.
2. Add the Julia [packages](http://docs.julialang.org/en/release-0.4/manual/packages/) Whorl depends on. The base code&#x2014;everything except `examples.jl`&#x2014;uses two packages:
  
  * [ValidatedNumerics](http://dpsanders.github.io/ValidatedNumerics.jl/)
  * [Compose](http://composejl.org/)
  
  If you want to run `examples.jl` as is, you'll need a few more packages:
  
  * [Colors](https://github.com/JuliaGraphics/Colors.jl)
  * [Gadfly](https://github.com/dcjones/Gadfly.jl)
  * [DataFrames](http://juliastats.github.io/DataFrames.jl/stable/)
  
  To [add](http://docs.julialang.org/en/release-0.4/manual/packages/#adding-and-removing-packages) packages, call the `Pkg.add` function from within Julia's interactive environment, which you can launch with the command-line call `julia`.

## Examples

You can use the examples in `examples.jl` to check that Whorl is working and to see what it can do. To run the examples from within Julia's [interactive environment](http://docs.julialang.org/en/release-0.4/manual/interacting-with-julia/), start with the following steps:

* Open a command line shell.
* Go to the folder you're keeping the Whorl source code in.
* Use the `julia` command to launch the interactive environment.
* Call `include("examples.jl")`. After some time, the environment should respond with the text `Examples`, showing that the `Examples` module is loaded.

Now you're ready to run the example functions described below. Keep in mind that Julia is just-in-time compiled, and it does a lot of caching too, so a function might take a few seconds to run the first time you call it.

### Abelianization

Call `Examples.abelianization_ex()`. The environment will print a description of an SL<sub>2</sub> **C** cocycle over an interval exchange and its abelianization. You can verify that the abelianized cocycle splits as a product of **C**<sup>&#x00d7;</sup> cocycles by observing that its transition matrices are diagonal.

### Shear parameter plots

### Geodesic lamination movie

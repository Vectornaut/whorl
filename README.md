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

Calling `Examples.abelianization_ex()` prints a description of an SL<sub>2</sub> **C** cocycle over an interval exchange and its abelianization. You can verify that the abelianized cocycle splits as a product of **C**<sup>&#x00d7;</sup> cocycles by observing that its transition matrices are diagonal.

The interval exchange used in this example is a first-return map for the vertical flow on a genus-5 translation surface *&Sigma;*, which happens to be the translation double cover of a genus-2 half-translation surface *C*. The cocycle over the interval exchange comes from the holonomy representation of a hyperbolic structure on *C*. The diagonal entries of the abelianized transition matrices are shear coordinates of the hyperbolic structure.

### Shear parameter plots

Calling `Examples.shear_plot_ex()` makes a plot showing how the shear coordinates described in the abelianization example depend on the vertical direction of *&Sigma;*. It also plots the generalized shear coordinates of two perturbed versions of the representation used in the abelianization example. The real and imaginary parts of the coordinates are plotted separately, in different colors. The coordinates of the unperturbed representation are real, as expected.

The plots appear in the working directory, in files called `no-perturbation.pdf`, `small-perturbation.pdf`, and `large-perturbation.pdf`. By default, `shear_plot_ex` produces very low-resolution plots, keeping the computation under 20 seconds on my laptop. You can get high-resolution plots by calling `Examples.shear_plot_ex(highres = true)`. The computation takes about 5 minutes on my laptop, but the results are worth it.

You can try out different perturbations by changing the array `perturbation` in the source code for `shear_plot_ex` and calling `include("examples.jl")` again. This reloads the module `Examples`, redefining the functions it provides.

### Geodesic lamination movie

The surface *C* described in the abelianization example comes with both a half-translation structure and a hyperbolic structure. The vertical foliation of the half-translation structure &ldquo;snaps tight&rdquo; to a geodesic lamination with respect to the hyperbolic structure. Calling `Examples.movie()` draws that geodesic lamination on the universal cover of *C*. More precisely, it draws the four ideal triangles that make up the complement of the lamination, with a different color for each triangle. The picture appears in the working directory, in the file `triangle_test.pdf`.

Calling `Examples.movie(testframe = false)` varies the vertical direction of *C* and draws the geodesic lamination for each direction. The pictures appear in the subdirectory `triangle_mov`. If you have [ImageMagick](http://www.imagemagick.org/script/convert.php) installed, running the script `splice` from inside `triangle_mov` splices the pictures together into an animated GIF.

### Looking under the hood

You can learn about the functions that make up Whorl by calling them from the interactive environment and examining the objects they return. For example, let's say you want to know more about the hyperbolic structure on *C* that features in all of the examples above. The holonomy representation of this hyperbolic structure enters our computations through a call to `Regular.generators(2)`, which returns a list of matrices. You can see the list by calling

    map(Examples.prettyprint, Regular.generators(2))

from the interactive environment. If the module prefixes get annoying to type, you can remove the need for them by making `using` statements from the interactive environment.

# 2d-induction-code
A simulation mostly written in FORTRAN which solves a 2d induction-like equation for a scalar field in a periodic Cartesian plane.

## The equation
The equation solved in this code is
<br>
<p align = "center">
<img src = "https://latex.codecogs.com/png.latex?\bg_white\frac{\partial%20T}{\partial%20t}%20=%20-\left(%20\frac{\partial}{\partial%20x}(v_xT)%20+%20\frac{\partial}{\partial%20y}(v_yT)%20\right)%20+%20\eta\left(%20\frac{\partial^2}{\partial%20x^2}%20+%20\frac{\partial^2}{\partial%20y^2}%20\right)T" alt = "equation1"/>
</p>
which is a special case of the MHD induction equation
<br>
<p align = 'center'>
<img src = "https://latex.codecogs.com/png.latex?\bg_white\frac{\partial\mathbf{B}}{\partial%20t}=\nabla\times(\mathbf{v}\times\mathbf{B}-\eta\nabla\times\mathbf{B})" alt = "equation2">
</p>
For more details, look into the document.

## The codes
This repository contains 3 program files:
1. `trid_solve.f90`: This program contains 3 subroutines: `trid`, `inner` and `cyclic_trid` (in this order).
   1. `trid`: solves a linear system of equations forming a tridiagonal matrix, using the Thomas algorithm.
   2. `inner`: calculates the inner product of two vectors.
   3. `cyclic_trid`: solves a linear system of equations forming a cyclic triadiagonal matrix, using the Sherman-Morrison formula and the Thomas algorithm.
   
   More details on the subroutines have been provided in the program file itself.

2. `induc_2d.f90`: This program solves the 2d induction-like equation using the subroutine `cyclic_trid` from `trid_solve.f90`. The algorithm has been explained in the document.
   One can modify the diffusivity, the velocity fields and the initial conditions, which are defined in the first half of the code. The default initial condition is a Gaussian.
   The code produces 4 document files, named `"induc_2d_dat_<num>"`, `"induc_2d_Xvals_<num>"`, `"induc_2d_Yvals_<num>"` and `"induc_2d_eta_t_<num>"`, where `<num>` is a number the user can put in the name before executing the code. These names are also defined at the beginning of the code as strings, named `filename_1` through `4` respectively.

3. `induc_2d_plot.py`: Python code to plot the outputs of `induc_2d.f90`. In the first part of the code, the user can change the names of the input (i.e., output of `induc_2d.f90`) and output (a `.gif` file) files of the code, and also the title of the plot.

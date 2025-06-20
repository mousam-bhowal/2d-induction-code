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
<img src = "https://latex.codecogs.com/png.latex?\bg_white\frac{\partial \mathbf{B}}{\partial t} = \nabla\times(\mathbf{v}\times\mathbf{B} - \eta\nabla\times\mathbf{B})" alt = "equation2">
</p>

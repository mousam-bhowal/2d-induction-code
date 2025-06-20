# 2d-induction-code
A simulation mostly written in FORTRAN which solves a 2d induction-like equation for a scalar field in a periodic Cartesian plane.

## The equation
The equation solved in this code is
$$
\frac{\partial T}{\partial t} = -\left( \frac{\partial}{\partial x}(v_xT) + \frac{\partial}{\partial y}(v_yT) \right) + \eta\left( \frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2} \right)\eta
$$

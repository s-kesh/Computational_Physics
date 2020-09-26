# To Solve ODE
Just define your f(x,y) (dy/dx = f(x,y) in function.f90
If it is higher order ODE, please define f(x,y) as a column vector

## Compile
```
	gfortran function.f90 ode.f90 diff.f90 -o <output>
```
For second part of assignment replace `diff.f90` to `shm.f90`

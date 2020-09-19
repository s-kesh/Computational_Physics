# To Solve ODE
Just define your f(x,y) (dy/dx = f(x,y) in diff.f90
If it is higher order ODE, please define f(x,y) as a column vector

## Compile
```
	gfortran ode.f90 diff.f90 -o <output>
```

## Running
Just run the compiled binary. It will ask user input itself.

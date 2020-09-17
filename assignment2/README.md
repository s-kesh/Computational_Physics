# To Solve ODE
Just define your f(x,y) (dy/dx = f(x,y) in function.f90
Edit diff.f90 to call correct method

## Compile
```
	gfortran function.f90 ode1d.f90 diff.f90 -o <output>
```

## Running
Just run the compiled binary. It will ask user input itself.

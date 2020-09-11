# To solve Particle in a Finite Box.
Just edit function and its diffentiation in source code.

## Compile
```
	gfortran <input> -o <output>
```

##Running

### For finding single eigenvalue in given segment
Just run the compiled binary. It will ask user input itself.

### For finding all possible eigenvalue of bound states
```
	bash runthis <output> <potential in eV> <half width in A⁰>
```

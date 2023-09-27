# SolitonSimulator

Non-linear optics involves the study of the interaction of light within non-linear media. The propagation of light through an optical fibre is described by Maxwell's equations and by considering the Polarization vector $\textbf{P}$ as a non-linear equation of the electric field $\textbf{E}$, which allows you to rewrite the Electromagnetic Wave equation to incorporate the non-linear behaviour. This can then be reduced to a form equivalent to the well known non-linear Schrodinger equation, the solution of which constitutes the focus of this project.

$$i \frac{\partial \Psi (z,t)}{\partial z} = -\frac{\partial^2 \Psi(z,t)}{\partial t^2} + c| \Psi(z,t)|^2 \Psi(z,t)$$

Numerical analysis using the split step method and RK4 to solve the non-linear Schrodinger equation were used to confirm the existence of the so called soliton solutions, which retain their shape and velocity during propagation, and are known to obey the non-linear Schrodinger equation. Using this scripit it is also possible to numerically confirm the existence of forces between neighbouring solitons that lead to their relative attraction in an oscillatory nature.

The python script contains two functions pertaining to the implementation of the RK4 method and the Split Step method for solving the NLSE. To switch between them, simpoly change the function on line 84.

# Example Output

![Figure_1](https://github.com/Matthew-Hill2000/SolitonSimulator/assets/124274792/cc50d2a0-69a2-4c9c-9980-94f50080146f)

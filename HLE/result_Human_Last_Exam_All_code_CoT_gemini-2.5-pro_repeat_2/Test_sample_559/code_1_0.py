# The separatrix of the system is a curve described by an equation relating d and u.
# Through analysis of the system's phase portrait, the separatrix is identified as an
# invariant manifold passing through a saddle point.
# This curve was found to be d = -u^2.

# We define the coefficient and the power for this equation.
coefficient = -1
power = 2

# Now, we print the final equation. The numbers that define the equation,
# the coefficient and the power, are explicitly part of the output string.
print(f"The equation for the separatrix is d = ({coefficient}) * u**({power})")
# The analysis of the system of differential equations reveals that a key separatrix
# is the unstable manifold of the saddle point at (d, u) = (-1, 1).
# We have found an analytical equation for this curve.

# The equation is of the form: d = a * u^p
# Let's define the numerical values for the parameters 'a' and 'p'.

# The coefficient of the u term
a = -1
# The exponent of the u term
p = 2

# Now, we print the final equation for the separatrix.
# The numbers that form the equation are explicitly printed as requested.
print("The equation of the separatrix is:")
print(f"d = {a} * u^{p}")
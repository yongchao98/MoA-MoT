# Based on the analysis of the system of differential equations,
# the separatrix is found by identifying an invariant manifold that
# passes through the saddle point of the system, which is located at (u, d) = (1, -1).
# The equation for this separatrix is d = -u^2.

# We define the components of this equation to print them clearly.
# The equation is of the form: d = coefficient * u^power
coefficient = -1
variable_name = "u"
power = 2

# We now print the equation for the separatrix, showing each number involved.
print(f"The equation for the separatrix is:")
print(f"d = {coefficient} * {variable_name}**{power}")
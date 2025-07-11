# In the path integral formalism for fermions, Grassmann variables are used.
# The properties of integration over these variables are defined to uphold
# the Pauli exclusion principle.

# The fundamental property of a Grassmann variable 'eta' (η) is that its square is zero: η² = 0.
# This directly models the Pauli principle, where a quantum state cannot be occupied
# by two identical fermions.

# The integration over these variables is defined by the Berezin integral.
# The "measure" of this integral doesn't have a value in itself, but it is defined
# by the value it yields when operating on a function. The key defining relation
# which normalizes the integral is:
#
# ∫ η dη = 1
#
# This code will print out this fundamental equation.

integral_value = 1

# Print the final equation with its components
print("The normalization of the measure for a Grassmann variable integral is defined by the following relation:")
print("Integral symbol: ∫")
print("Measure: dη")
print("Integrand: η")
print("Result: {}".format(integral_value))
print("\nWhich forms the equation:")
print("∫ η dη = {}".format(integral_value))
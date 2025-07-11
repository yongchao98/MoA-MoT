# The Pauli exclusion principle is enforced in fermionic path integrals
# through the use of anticommuting Grassmann variables. The fundamental
# integral, which defines the measure, is the Gaussian Grassmann integral.
# This script prints the equation for this integral and its result.

# Define the symbolic parts of the equation
integral = "∫dη̄ dη"
exponent = "e^(-η̄aη)"
equals = "="
result = "a"

# Print the full equation.
# The value 'a' is the result of the integral, representing the measure.
print("The equation for the fundamental Grassmann integral is:")
print(f"{integral} {exponent} {equals} {result}")
print("\nWhere:")
print(" - η, η̄ are anticommuting Grassmann variables whose properties encode the Pauli exclusion principle.")
print(f" - The value of the integral is '{result}', the coefficient from the exponent.")

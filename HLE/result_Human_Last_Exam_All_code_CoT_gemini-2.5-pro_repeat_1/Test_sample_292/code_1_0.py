import math

# The size of the vocabulary and the length of the sequences.
n = 99

# Based on the mathematical derivation, the sum S is equal to (n+1)^(n-1).
# First, let's calculate the base and exponent of this expression.
base1 = n + 1
exponent1 = n - 1

# The problem asks for the answer to be written as a power of 10.
# We can rewrite the expression: (10^2)^98 = 10^(2*98) = 10^196.
base2 = 10
exponent2 = 2 * exponent1

# The final equation is 100^98 = 10^196.
# As requested, we will print each number in this final equation.
print(f"The calculation is based on the simplified formula: ({n}+1)^({n}-1)")
print(f"This evaluates to: {base1}^{exponent1}")
print(f"Expressed as a power of 10, this is: {base2}^{exponent2}")
print("\nThe numbers in the final equation are:")
print(base1)
print(exponent1)
print(base2)
print(exponent2)

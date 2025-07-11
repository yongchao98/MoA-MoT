import math

# Based on the derivation, the maximal overhang is (2 + sqrt(2))/3 block sidelengths.
# We match this to the format (a + sqrt(b))/(1+c).

# Numerator: a + sqrt(b) = 2 + sqrt(2)
a = 2
b = 2

# Denominator: 1 + c = 3
c = 2

# Print the determined integers.
print(f"The solution is a={a}, b={b}, c={c}.")

# To verify, calculate the value of the expression
overhang = (a + math.sqrt(b)) / (1 + c)
print(f"The maximal overhang is ({a} + sqrt({b}))/(1 + {c}) = {overhang:.4f} block sidelengths.")

# The final answer format requires printing the space-separated integers.
print("\nThe final integers a, b, c are:")
print(a, b, c)
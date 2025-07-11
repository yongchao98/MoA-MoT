import sys

# This script is for demonstration purposes to output the result symbolically.
# Ordinal arithmetic is not a standard feature in Python.

# The alphabet {a, b, c, d} has an order type of 4.
d = 4

# The order type of finite strings over an alphabet with order type d
# under lexicographical ordering is omega^d.
# We will represent omega symbolically with the character 'ω'.

# The result is ω^4. We can represent this as a polynomial in ω.
# ω^4 = 1*ω^4 + 0*ω^3 + 0*ω^2 + 0*ω + 0
print("The order type is omega^d where d is the order type of the alphabet.")
print("The alphabet {a,b,c,d} has order type d = 4.")
print("So the final order type is omega^4.")
print("\nThis can be represented as an ordinal polynomial:")
# Printing the 'equation' part as requested
print("1 * ω^4 + 0 * ω^3 + 0 * ω^2 + 0 * ω^1 + 0")

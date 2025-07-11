# The graph G is a specific type of decorated ladder graph.
# The connective constant (μ) of an infinite graph is the asymptotic growth rate
# of the number of self-avoiding walks. For this particular graph, the
# connective constant is known to be μ = sqrt(2 + sqrt(5)).

# To find the minimal polynomial of μ over the rational numbers, we perform
# algebraic manipulation. Let x = μ.
# x = sqrt(2 + sqrt(5))
# Square both sides:
# x^2 = 2 + sqrt(5)
# Isolate the square root term:
# x^2 - 2 = sqrt(5)
# Square both sides again to eliminate the root:
# (x^2 - 2)^2 = 5
# Expand the equation:
# x^4 - 4x^2 + 4 = 5
# This simplifies to the polynomial equation:
# x^4 - 4x^2 - 1 = 0

# This polynomial, P(x) = x^4 - 4x^2 - 1, is the minimal polynomial because
# it is irreducible over the rational numbers.

# We represent the final polynomial equation x^4 + 0*x^3 - 4*x^2 + 0*x - 1 = 0
# by its integer coefficients.
# The following code prints each of these coefficients as requested.

c4 = 1  # Coefficient of x^4
c3 = 0  # Coefficient of x^3
c2 = -4 # Coefficient of x^2
c1 = 0  # Coefficient of x^1
c0 = -1 # Coefficient of x^0 (the constant term)

print("The numbers in the final equation P(x)=0 are the coefficients of the polynomial:")
print(c4)
print(c3)
print(c2)
print(c1)
print(c0)

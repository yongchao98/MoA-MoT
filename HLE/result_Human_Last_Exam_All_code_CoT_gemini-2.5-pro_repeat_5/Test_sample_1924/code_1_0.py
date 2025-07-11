import sympy

# Define the symbols for the lengths of the segments
MA, MB, MG, MH = sympy.symbols('MA MB MG MH')

# The problem asks to express MG - MH in terms of MA and MB.
# Based on geometric theorems related to the power of a point and concyclic points
# (specifically, a generalization of the Butterfly Theorem or Haruki's Theorem),
# it can be proven that the following relationship holds:
# vec(MA) + vec(MB) = vec(MG) + vec(MH)
# where vec represents signed distances from M.

# Let M be the origin on the line AB. Let the direction from M to B be positive.
# The coordinate of A is -MA.
# The coordinate of B is +MB.
# Let the coordinates of G and H be g and h.
# The theorem states: (-MA) + (+MB) = g + h.

# We need to find the value of MG - MH (difference in lengths).
# It can be shown that one of G or H lies on the segment MA and the other on MB.
# Let's assume G is on the side of A (g < 0) and H is on the side of B (h > 0).
# Then the length MG = |g| = -g.
# And the length MH = |h| = h.
# So, MG - MH = -g - h = -(g + h).

# From the theorem, g + h = MB - MA.
# Therefore, MG - MH = -(MB - MA) = MA - MB.

# We will now create and print the equation representing this result.
# The coefficients are 1 for MA and -1 for MB.
equation = sympy.Eq(MG - MH, 1*MA - 1*MB)

# Print the final equation
print("The relationship is:")
print(f"{equation.lhs} = {equation.rhs}")

print("\nExplaining the coefficients in the final equation:")
print(f"The coefficient of MA is: {1}")
print(f"The coefficient of MB is: {-1}")

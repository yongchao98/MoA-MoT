# This script provides the solution to the mathematical question about the
# property (Rn) of Noetherian local rings.

# The question asks for the smallest non-negative integer n such that the property
# (Rn), meaning "regular in codimension n", is not preserved by completion.

# Based on established theorems in commutative algebra:
# - For n=0, the property (R0) is preserved.
# - For n=1, the property (R1) is preserved.
# - For n=2, the property (R2) is NOT preserved, as shown by a
#   counterexample constructed by Masayoshi Nagata.

# Therefore, the smallest such integer is 2.
# We can state this result as the simple equation:
# n = 2

# The final equation is n = 2. We print the number in this equation.
final_answer = 2
print("The smallest non-negative integer n such that (Rn) is not preserved by completion is:")
print(final_answer)
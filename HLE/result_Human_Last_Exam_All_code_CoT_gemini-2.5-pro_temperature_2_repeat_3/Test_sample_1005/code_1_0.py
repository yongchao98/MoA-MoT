import math

# The problem simplifies to calculating 2^(15/16).
# The final equation is result = A^(B/C)
A = 2
B = 15
C = 16

# As per the instruction, we output each number in the final equation.
print("A =", A)
print("B =", B)
print("C =", C)

result = A ** (B / C)
print(f"The result of the definite integral is A**(B/C) = {A}**({B}/{C})")
print(result)

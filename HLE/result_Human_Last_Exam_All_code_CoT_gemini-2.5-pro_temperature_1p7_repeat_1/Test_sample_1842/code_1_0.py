# The problem reduces to solving x^4 + y^4 + z^4 = w^4.
# We are looking for the solution (x, y, z, w) with the smallest
# maximum value among its components. This is the smallest known
# counterexample to Euler's sum of powers conjecture for n=4.

# Assign the values from this known solution.
# The variables x, y, z can be any permutation of these three numbers.
x = 95800
y = 217519
z = 414560
# w would be 422481, but it's not needed for the sum x+y+z.

# Calculate the required sum
final_sum = x + y + z

# Print the calculation as requested.
# The 'f' string formatting allows embedding variables directly into the string.
print(f"The equation to solve is {x}^4 + {y}^4 + {z}^4 = 422481^4.")
print(f"The sum x + y + z is calculated as:")
print(f"{x} + {y} + {z} = {final_sum}")
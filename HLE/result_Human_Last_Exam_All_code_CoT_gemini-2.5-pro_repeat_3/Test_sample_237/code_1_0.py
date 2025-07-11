import math

# Step 1: Determine the value of z.
# Based on the analysis of XVM data types and the algorithm to calculate 100!,
# the smallest memory size 'z' is calculated as follows:
# - result array: An array of 158 'digit' types = 158 * 1D = 158D
# - Loop counter 'i' (2 to 100): 'char' type = 3D
# - Inner loop counter 'j' (0 to 157): 'char' type = 3D
# - Number of digits tracker 'num_digits': 'char' type = 3D
# - Carry variable for multiplication: 'char' type = 3D
# Total z = 158 + 3 + 3 + 3 + 3 = 170.
z = 170

# Step 2: Calculate the value of y.
# Use Python's math library to compute 100!
factorial_100 = math.factorial(100)

# Convert the result to a string and get the first 3 characters.
y_string = str(factorial_100)[:3]
y = int(y_string)

# Step 3: Output the result in the format z:y.
# The instruction "output each number in the final equation!" means we should
# print the components of the final answer 'z:y'.
print(z, end="")
print(":", end="")
print(y)
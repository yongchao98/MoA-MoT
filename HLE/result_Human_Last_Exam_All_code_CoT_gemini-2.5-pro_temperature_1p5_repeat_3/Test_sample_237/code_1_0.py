import math

# This script calculates the necessary values and prints the final answer
# in the required z:y format.

# Step 1: Calculate 'z', the minimal memory size in decimal digits (D).
# Based on the analysis of the most memory-efficient C program for the Wuxing VM,
# the sizes of the required variables are:
res_size = 158  # 158 * 1D for the result array
i_size = 3      # 1 * 3D for the char loop counter 'i'
j_size = 3      # 1 * 3D for the char loop counter 'j'
size_tracker_size = 3 # 1 * 3D for the char size tracker
carry_size = 2  # 1 * 2D for the cent carry variable
product_size = 3# 1 * 3D for the char temporary product variable

# The total memory size z is the sum of these components.
z = res_size + i_size + j_size + size_tracker_size + carry_size + product_size

# Step 2: Calculate 'y', the first 3 digits of 100!
# Python's math.factorial can handle large numbers, making this straightforward.
factorial_100 = math.factorial(100)
factorial_100_str = str(factorial_100)
y = factorial_100_str[0:3]

# Step 3: Print the final answer in the format z:y
# The prompt asks to output each number in the final equation.
# Here, the equation is the z:y format, and we output the calculated z and y.
print(f"{z}:{y}")

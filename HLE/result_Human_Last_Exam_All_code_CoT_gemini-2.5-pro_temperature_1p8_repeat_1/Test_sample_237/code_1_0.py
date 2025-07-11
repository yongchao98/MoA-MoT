import math

# Part 1: Calculate 'z', the smallest total memory size in decimal digits (D).
# This is based on designing an optimized C program for the Wuxing architecture.
# Data type sizes in Digits (D): digit=1, char=3.

# Memory for the array to hold the 158 digits of 100!.
# We use an array of 160 'digit's for safety and efficiency on the decimal machine.
mem_result_array = 160 * 1 # 160 D

# Memory for the main loop counter 'i' (from 2 to 100).
# The 'char' type (range 0-999) is the smallest suitable type.
mem_i = 3 # 3 D

# Memory for the inner loop counter 'j' (iterating up to 159).
# 'char' is also the smallest suitable type here.
mem_j = 3 # 3 D

# Memory for tracking the number of digits in the result (up to 158).
# 'char' is the smallest suitable type.
mem_num_digits = 3 # 3 D

# Memory for the 'carry' variable during multiplication.
# The max value is ~990 (9 * 100 + 90), which fits within a 'char'.
mem_carry = 3 # 3 D

# 'z' is the sum of memory for all variables.
z = mem_result_array + mem_i + mem_j + mem_num_digits + mem_carry

# Part 2: Calculate 'y', the first 3 digits of 100!.
# Use Python's math library for an exact value.
factorial_100 = math.factorial(100)

# Convert to string and slice the first three characters.
y = str(factorial_100)[:3]

# Part 3: Print the final output.
# As requested, we first output the numbers in the final equation for 'z'.
print(f"{z} = {mem_result_array} + {mem_i} + {mem_j} + {mem_num_digits} + {mem_carry}")

# Then, print the final answer in the format z:y.
print(f"{z}:{y}")
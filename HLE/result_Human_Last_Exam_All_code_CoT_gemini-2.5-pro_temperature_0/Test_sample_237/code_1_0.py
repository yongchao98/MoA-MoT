import math

# This script calculates the answer for the Wuxing VM problem.
# The answer is in the format z:y, where:
# z = smallest memory size in decimal digits (D) for an optimized C program.
# y = the first 3 digits of 100!

# --- Step 1: Calculate the smallest memory size 'z' ---

# We design an optimized C program using BigInt arithmetic with a base of 100.
# This requires an array to hold the digits of the large number.
#
# XVM Data type sizes in decimal digits (D):
# cent: 2D
# int: 6D
# unsigned int: 6D
#
# The variables needed for the C program and their memory cost in D are:
# - result[79]: An array of 'cent' type to store the 158-digit result (79 * 2D).
# - i: An 'int' loop counter for the factorial (2 to 100). 'cent' (max 99) is too small.
# - j: A 'cent' loop counter for the array (up to 78).
# - size: A 'cent' to track the array's used size (up to 79).
# - carry: A 'cent' for multiplication carry-over (max carry is 99).
# - product: An 'unsigned int' temporary variable for intermediate products
#   (e.g., 99 * 100 + 99 = 9999, which fits in an unsigned int).

mem_result_array = 79 * 2
mem_i_counter = 6
mem_j_counter = 2
mem_size_tracker = 2
mem_carry = 2
mem_product = 6

# 'z' is the sum of the memory used by all variables.
z = mem_result_array + mem_i_counter + mem_j_counter + mem_size_tracker + mem_carry + mem_product

# --- Step 2: Calculate the first 3 digits 'y' of 100! ---

factorial_100 = math.factorial(100)
y = str(factorial_100)[:3]

# --- Step 3: Output the result ---

# The prompt requires outputting each number in the final equation for 'z'.
print("Calculation for z (total memory in D):")
print(f"z = {mem_result_array} (result[79]) + {mem_i_counter} (i) + {mem_j_counter} (j) + {mem_size_tracker} (size) + {mem_carry} (carry) + {mem_product} (product)")

# Print the final answer in the format z:y
print("\nFinal Answer:")
print(f"{z}:{y}")
# This script calculates and prints the final answer based on the analysis.
# z: The smallest memory size in decimal digits (D) for an optimized C program
#    to calculate 100! on the XVM architecture.
# y: The first 3 digits of the result of 100!.

# Calculation for z (memory size in D)
# The result of 100! (158 digits) is stored in an array of 79 'cent' types.
size_result_array = 79 * 2  # 158D
# The loop counter 'i' (up to 100) requires a 'char' type.
size_i = 3
# The inner loop counter 'j' (up to 78) fits in a 'cent' type.
size_j = 2
# A variable to track the array size (up to 79) fits in a 'cent' type.
size_result_size = 2
# The 'carry' variable in multiplication will not exceed 99, so a 'cent' is sufficient.
size_carry = 2
# The intermediate 'product' can be up to 9999, requiring an 'int' type.
size_product = 6

z = size_result_array + size_i + size_j + size_result_size + size_carry + size_product

# Value for y (first three digits of 100!)
# 100! starts with 933...
y = 933

# Print the final answer in the format z:y
print(f"{z}:{y}")
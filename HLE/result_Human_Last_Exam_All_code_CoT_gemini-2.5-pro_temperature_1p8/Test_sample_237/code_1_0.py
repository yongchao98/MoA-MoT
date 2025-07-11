import math

# Step 1: Determine 'z', the smallest memory size in D.
# This value is derived from designing an optimized C program for the hypothetical XVM.
#
# - Bignum storage: The result of 100! has 158 digits.
#   Using an array of 'unsigned char' (3D, base 1000) is most memory-efficient.
#   Array size: ceil(158 / 3) = 53 elements.
#   Memory for array: 53 * 3D = 159D.
#
# - Loop variables:
#   - Outer loop (i from 2 to 100): needs 'unsigned char' (3D).
#   - Inner loop (j from 0 to 52): needs 'unsigned cent' (2D).
#
# - Calculation variable:
#   - Intermediate product (`d[j] * i + carry`) can reach ~99,999.
#     This requires an 'unsigned int' (6D).
#
# - Total Memory (z):
#   z = 159 (array) + 3 (i) + 2 (j) + 6 (product) = 170D.
z = 170

# Step 2: Determine 'y', the first 3 digits of 100!.
# Python's math.factorial can compute this large number directly.
factorial_100 = math.factorial(100)

# Convert the integer result to a string to easily access the first digits.
y = str(factorial_100)[:3]

# Step 3: Print the final answer in the required format 'z:y'.
# As per the instruction, the code prints the value of each number in the final result.
print(f"{z}:{y}")
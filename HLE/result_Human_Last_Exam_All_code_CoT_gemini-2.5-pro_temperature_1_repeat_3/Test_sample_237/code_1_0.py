import math

# Step 1: Determine the number of digits in 100! to size the result array.
# The result of 100! must be stored. We'll use an array of digits for this.
factorial_100 = math.factorial(100)
# The number of digits determines the minimum size of our result array.
num_digits_in_result = len(str(factorial_100))

# The most memory-efficient type to store a single digit is 'digit' (1D).
size_result_array = num_digits_in_result * 1  # Size in D

# Step 2: Determine the minimum data type size for all other necessary variables
# in a C program that calculates the factorial.

# Variable 'i' (for the main loop: for i = 2 to 100):
# Needs to hold 100. 'cent' (0-99, 2D) is too small.
# The next smallest type is 'char' (0-999, 3D).
size_i = 3  # D

# Variable 'j' (for the inner loop, iterating through the digits array):
# The index will go up to num_digits_in_result - 1 (157).
# 'char' (0-999, 3D) is the smallest type that can hold 157.
size_j = 3  # D

# Variable 'carry' (for multiplication):
# In an operation like (digit * 100 + carry), the max carry is 99.
# E.g., (9 * 100 + 99) / 10 = 99.9. Max carry is 99.
# 'cent' (0-99, 2D) is a perfect fit.
size_carry = 2  # D

# Variable 'num_digits' (to track the current size of the result):
# This will count up to 158.
# 'char' (0-999, 3D) is the smallest type that can hold 158.
size_num_digits = 3 # D

# Temporary variable 'temp_prod' (to hold intermediate products):
# Stores the result of `result[j] * i + carry`.
# Maximum value: 9 * 100 + 99 = 999.
# 'char' (0-999, 3D) is a perfect fit.
size_temp_prod = 3 # D

# Step 3: Sum the sizes to get z, the total minimum memory usage.
z = size_result_array + size_i + size_j + size_carry + size_num_digits + size_temp_prod

# Step 4: Get y, the first 3 digits of 100!.
y = str(factorial_100)[:3]

# Step 5: Print the final answer in the format z:y.
print(f"{z}:{y}")
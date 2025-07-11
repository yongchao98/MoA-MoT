import math

# Plan to calculate z, the smallest memory size in D for an optimized C program
# that computes 100! on the Wuxing Virtual Machine (XVM).

# 1. Determine the number of digits in 100! to size the result array.
# The result of 100! is a 158-digit number.
# The most memory-efficient way to store this is an array of 'digit' (1D) type.
num_digits_in_result = 158
memory_for_result_array = num_digits_in_result * 1  # 1D for each 'digit'

# 2. Select the smallest data type for each required variable in the factorial algorithm.
# Data Types: digit(1D, 0-9), cent(2D, 0-99), char(3D, 0-999)

# - Loop counter 'i' from 2 to 100.
#   A 'cent' (0-99) is too small. 'char' (0-999) is the next smallest.
memory_for_i = 3  # D

# - Inner loop counter 'j' for digits, up to 157.
#   A 'cent' (0-99) is too small. 'char' (0-999) is required.
memory_for_j = 3  # D

# - Variable 'size' to track the number of digits, up to 158.
#   A 'cent' (0-99) is too small. 'char' (0-999) is required.
memory_for_size = 3 # D

# - 'carry' variable for multiplication. Max value of carry is 99.
#   A 'cent' (0-99) is a perfect fit.
memory_for_carry = 2  # D

# - 'temp' variable for holding the intermediate product (digit * i + carry).
#   Max value is 9 * 100 + 99 = 999. 'char' (0-999) is a perfect fit.
memory_for_temp = 3  # D

# 3. Sum the memory of all variables to get z.
z = memory_for_result_array + memory_for_i + memory_for_j + memory_for_size + memory_for_carry + memory_for_temp

# 4. Calculate y, the first 3 digits of 100!.
factorial_100 = math.factorial(100)
y = int(str(factorial_100)[:3])

# 5. Print the final answer in the format z:y.
print(f"{z}:{y}")

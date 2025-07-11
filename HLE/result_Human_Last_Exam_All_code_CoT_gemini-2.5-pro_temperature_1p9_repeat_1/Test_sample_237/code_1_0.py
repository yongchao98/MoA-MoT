import math

# Step 1: Calculate 'y', the first three digits of 100!.
# This involves computing the factorial, converting it to a string, and slicing.
factorial_100 = math.factorial(100)
y = str(factorial_100)[:3]

# Step 2: Calculate 'z', the minimum memory size in decimal digits (D).
# This is based on the variables required for an optimized C program on XVM.
# 100! has 158 digits, so we need an array for arbitrary-precision arithmetic.
num_digits_in_result = len(str(factorial_100)) # This evaluates to 158

# Analysis of variables for the algorithm:
# - A `result` array to store the digits of 100!.
#   - Size: Needs to hold 158 digits.
#   - Type: An array of `digit` type. A `digit` is 1D.
#   - Memory cost: 158 * 1D = 158D.
mem_result_array = num_digits_in_result * 1

# - A loop counter `i` from 2 to 100.
#   - Max value: 100.
#   - Type: `cent` (0-99) is too small. `char` (0-999) is the smallest fit.
#   - Memory cost: 3D.
mem_i = 3

# - A variable `size` to track the number of digits in the array.
#   - Max value: 158.
#   - Type: `char` (0-999) is the smallest fit.
#   - Memory cost: 3D.
mem_size = 3

# - An inner loop counter `j` to iterate through the digits.
#   - Max value: 157 (as an index).
#   - Type: `char` (0-999) is the smallest fit.
#   - Memory cost: 3D.
mem_j = 3

# - A variable `carry` for multiplication.
#   - Max value is `(9 * 100 + 99) / 10 = 99`.
#   - Type: `cent` (0-99) fits perfectly.
#   - Memory cost: 2D.
mem_carry = 2

# - A temporary variable `temp` to hold `(result[j] * i + carry)`.
#   - Max value is `9 * 100 + 99 = 999`.
#   - Type: `char` (0-999) fits perfectly.
#   - Memory cost: 3D.
mem_temp = 3

# Summing the memory costs gives 'z'.
z = mem_result_array + mem_i + mem_size + mem_j + mem_carry + mem_temp

# Step 3: Print the final answer in the format z:y.
print(f"{z}:{y}")
<<<172:933>>>
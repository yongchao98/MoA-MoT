import math

# Step 1: Calculate 100! to get its properties.
factorial_100 = math.factorial(100)
factorial_str = str(factorial_100)

# Step 2: Determine 'y', the first 3 digits of the result.
# The full value of 100! is 9332621544394415268...
y = factorial_str[0:3]

# Step 3: Determine 'z', the smallest memory size in D.
# We need to design an optimized C program and calculate its variable memory usage.

# Sizes of the Wuxing data types in decimal digits (D)
size_d = {
    'digit': 1,
    'cent': 2,
    'char': 3,
    'int': 6,
    'long': 12
}

# The result of 100! has 158 digits. We need an array to store it.
# The most memory-efficient way is an array of 'digit's.
num_digits_in_result = len(factorial_str)
mem_result_array = num_digits_in_result * size_d['digit']

# Memory for temporary variables needed for the calculation loop.
# We choose the smallest possible data type for each variable.
# for (i = 2; i <= 100; i++)
mem_i = size_d['char']  # max value 100, fits in char (0-999)

# for (j = 0; j < num_digits; j++)
mem_j = size_d['char']  # max value 157, fits in char (0-999)

# A variable to track the number of digits in the result array
mem_num_digits_tracker = size_d['char'] # max value 158, fits in char (0-999)

# In each multiplication step: temp = result[j] * i + carry;
# Max temp = 9 * 100 + 99 = 999.
mem_temp_product = size_d['char'] # max value 999, fits in char (0-999)

# carry = temp / 10;
# Max carry = 999 / 10 = 99.
mem_carry = size_d['cent'] # max value 99, fits in cent (0-99)

# Calculate z, the total memory for all variables.
z = mem_result_array + mem_i + mem_j + mem_num_digits_tracker + mem_temp_product + mem_carry

# Print the final result in the format z:y
# The numbers in the final equation are z and y.
print(f"{z}:{y}")
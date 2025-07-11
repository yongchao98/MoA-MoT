import math

# Part 1: Calculate 'z', the smallest memory size in D.
# This requires a theoretical analysis of the optimal C program for the XVM.
# The C program must calculate 100! using big number arithmetic, as 100! is a 158-digit number.
# We determine the memory for each variable based on the smallest possible data type.

# Variable 1: The array to store the big number.
# To store 158 digits, the most memory-efficient method is an array of `digit`s.
# We'll allocate a safe size of 160.
# In C: `digit result[160];`
# Memory for `result`: 160 elements * 1D/element = 160D.
mem_result = 160

# Variable 2: A counter for the number of digits currently in use.
# This counter will go up to 158. The `char` type (3D, range 0-999) is sufficient.
# In C: `char num_digits;`
# Memory for `num_digits`: 1 * 3D = 3D.
mem_num_digits = 3

# Variable 3: The main loop counter 'i' (from 2 to 100).
# The variable must hold the value 100. A `cent` (2D, range 0-99) is too small.
# The `char` type (3D, range 0-999) is the next smallest sufficient type.
# In C: `char i;`
# Memory for `i`: 1 * 3D = 3D.
mem_i = 3

# Variable 4: The inner loop counter 'j' (iterating through digits).
# This counter will go up to `num_digits - 1`, which is max 157.
# The `char` type (3D, range 0-999) is sufficient.
# In C: `char j;`
# Memory for `j`: 1 * 3D = 3D.
mem_j = 3

# Variable 5: The 'carry' value during multiplication.
# The max carry can be calculated as floor((9 * 100 + 99) / 10) = 99.
# The `cent` type (2D, range 0-99) is exactly sufficient.
# In C: `cent carry;`
# Memory for `carry`: 1 * 2D = 2D.
mem_carry = 2

# Variable 6: The temporary 'product' (result[j] * i + carry).
# The max product is 9 * 100 + 99 = 999.
# The `char` type (3D, range 0-999) is exactly sufficient.
# In C: `char product;`
# Memory for `product`: 1 * 3D = 3D.
mem_product = 3

# The total minimum memory 'z' is the sum of memory for all variables.
z = mem_result + mem_num_digits + mem_i + mem_j + mem_carry + mem_product

# Part 2: Calculate 'y', the first 3 digits of 100!
# Python's math.factorial can handle large numbers.
factorial_100 = math.factorial(100)
# Convert to string and slice the first 3 characters.
y = str(factorial_100)[:3]

# Print the final result in the format z:y.
# The two numbers in the final output equation are z and y.
print(f"{z}:{y}")
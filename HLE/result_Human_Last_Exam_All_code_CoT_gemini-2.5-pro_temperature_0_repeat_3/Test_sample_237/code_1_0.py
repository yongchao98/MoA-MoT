import math

# Step 1 & 2: Determine the memory size 'z' for the hypothetical C program.
# The result of 100! has 158 digits. We need an array to store these digits.
# We'll use an array of 'digit' types, which are 1D each. A size of 160 is safe.
# memory_for_array = 160 * 1D = 160D

# Variables needed for the calculation:
# - result_array[160]: An array of 160 'digit' types. Size = 160 * 1D = 160D.
# - i: Loop counter from 2 to 100. Needs to hold 100. Smallest type is 'char' (3D).
# - num_digits: Tracks number of digits in result. Max 158. Smallest type is 'char' (3D).
# - j: Inner loop counter for the array. Max 159. Smallest type is 'char' (3D).
# - product: Temp var for (digit * i + carry). Max value is ~9*100+carry < 1000. Smallest type is 'char' (3D).
# - carry: Temp var for carry-over. Max value is ~product/10 < 100. Smallest type is 'cent' (2D).

z_memory_size_d = 160 + 3 + 3 + 3 + 3 + 2

# Step 3: Calculate the first 3 digits 'y' of 100!
factorial_100 = math.factorial(100)
y_first_three_digits = int(str(factorial_100)[:3])

# Step 4: Print the final answer in the format z:y
# The final equation is the formatted string z:y
print(f"{z_memory_size_d}:{y_first_three_digits}")

# The final answer is 174:933
# z = 174
# y = 933
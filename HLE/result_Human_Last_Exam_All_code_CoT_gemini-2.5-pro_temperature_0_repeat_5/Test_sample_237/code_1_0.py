import math

# This script calculates the answer to the Wuxing 100! problem.
# It determines 'z' (minimum memory in D) and 'y' (first 3 digits of 100!).

# Step 1: Define the sizes of the Wuxing data types in decimal digits (D).
data_type_sizes = {
    'digit': 1,
    'cent': 2,
    'char': 3,
    'int': 6,
    'long': 12
}

# Step 2: Analyze the memory requirements for an optimized C program.
# The program must calculate 100!, a 158-digit number. This requires an array.

# Variable 1: An array to store the 158 digits of the result.
# Each element stores one digit, so the 'digit' type is optimal.
num_digits_in_result = len(str(math.factorial(100)))
mem_result_array = num_digits_in_result * data_type_sizes['digit']

# Variable 2: A loop counter 'i' for the factorial (from 2 to 100).
# Max value is 100. The smallest type that can hold 100 is 'char' (0-999).
mem_loop_i = data_type_sizes['char']

# Variable 3: A variable 'result_size' to track the current number of digits.
# Max value will be 158. The smallest type is 'char' (0-999).
mem_result_size = data_type_sizes['char']

# Variable 4: An inner loop counter 'j' for the multiplication function.
# Max value will be the array index, 157. The smallest type is 'char' (0-999).
mem_loop_j = data_type_sizes['char']

# Variable 5: A 'carry' variable for the multiplication logic.
# Analysis shows the max carry is around 189 (from floor((9*100 + 999)/10)).
# The smallest type that can hold 189 is 'char' (0-999).
mem_carry = data_type_sizes['char']

# Variable 6: An intermediate 'product' variable in the multiplication.
# Analysis shows the max product is around 1089 (from 9*100 + 189).
# 'char' is too small. The smallest type is 'int' (up to 999,999).
mem_product = data_type_sizes['int']

# Step 3: Calculate 'z' by summing the memory sizes.
z = mem_result_array + mem_loop_i + mem_result_size + mem_loop_j + mem_carry + mem_product

# Step 4: Calculate 'y' by finding the first 3 digits of 100!.
factorial_100_str = str(math.factorial(100))
y = factorial_100_str[:3]

# Step 5: Print the final answer in the required z:y format.
# The prompt asks to "output each number in the final equation!".
# This print statement shows the final numbers for z and y.
print(f"{z}:{y}")
import math

# Step 1: Determine the smallest memory size 'z' based on the plan.
# As calculated, the minimum memory required for an optimized program is 170D.
# z = size_of(result_array[158]) + size_of_3_char_vars + size_of_carry_char
# z = 158*1D + 3*3D + 1*3D
z = 158 + 9 + 3

# Step 2: Calculate 100! and find the first 3 digits for 'y'.
factorial_100 = math.factorial(100)
# Convert the number to a string to slice the first three characters.
y_str = str(factorial_100)[:3]
y = int(y_str)

# Step 3: Print the final answer in the format z:y.
# The prompt asks to "output each number in the final equation!".
# We will print the values of z and y that form the final answer.
print(f"{z}:{y}")
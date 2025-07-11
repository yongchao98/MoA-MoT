# Define the dimensions of the rectangle.
n = 1000000001
m = 1000000001

# Step 1: Calculate the memory usage for an optimal C program.
# The value of n*m is ~10^18, which exceeds the capacity of a 32-bit integer.
# In C, we would need to use 'long long int' to store n, m, and the result,
# as it's guaranteed to be at least 8 bytes and can hold values up to ~9*10^18.
# The problem asks for the memory used by variables for n, m, and the output.
size_of_long_long_int = 8  # in bytes
memory_for_n = size_of_long_long_int
memory_for_m = size_of_long_long_int
memory_for_output = size_of_long_long_int

# The total memory usage 'x' is the sum for these three variables.
x = memory_for_n + memory_for_m + memory_for_output

# Step 2: Calculate the maximal covered area.
# Each brick has an area of 2. The total covered area must be a multiple of 2.
# The maximal area is the largest even number less than or equal to n * m.
# This can be calculated with integer arithmetic.
# Python's '//' operator performs floor division.
total_area = n * m
max_covered_area = (total_area // 2) * 2

# Step 3: Output the result in the required format 'x:o'.
# The final equation is the one that forms the output string.
# We are printing the two numbers, x and max_covered_area, that make up this final "equation".
print(f"{x}:{max_covered_area}")

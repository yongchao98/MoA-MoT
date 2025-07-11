# Define the dimensions of the rectangle
n = 1000000001
m = 1000000001

# 1. Analyze memory usage for an equivalent C program.
# The values for N, M, and the result (N*M) require a 64-bit integer type
# (long long in C) to prevent overflow. A long long is 8 bytes.
# We need variables for n, m, and the output.
size_of_long_long_in_bytes = 8
memory_usage_x = size_of_long_long_in_bytes * 3

# 2. Calculate the maximal covered area.
# Since Python integers handle arbitrary size, we don't need to worry about overflow here.
total_area = n * m

# The maximal covered area is the largest even number less than or equal to the total area.
# This can be calculated as total_area - (total_area % 2).
maximal_area_o = total_area - (total_area % 2)

# 3. Print the final result in the format "x:o".
# The instruction "output each number in the final equation!" is interpreted as
# printing the two components of the final answer, 'x' and 'o'.
print(f"{memory_usage_x}:{maximal_area_o}")
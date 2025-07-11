# Define the dimensions of the rectangle as given in the problem.
N = 1000000001
M = 1000000001

# Calculate the maximal covered area.
# The area is the largest even number less than or equal to the total area N * M.
# This can be calculated using integer division: (N * M // 2) * 2.
# Python's arbitrary-precision integers handle the large numbers automatically.
output_area = (N * M // 2) * 2

# Calculate the memory usage for an equivalent optimal C program.
# The variables n, m, and the output would be stored.
# n = 1000000001 requires a 64-bit integer (long long in C).
# m = 1000000001 requires a 64-bit integer (long long in C).
# The product N*M requires a 64-bit unsigned integer (unsigned long long in C).
# On most systems, the size of these types is 8 bytes.
mem_n_bytes = 8
mem_m_bytes = 8
mem_output_bytes = 8

# The final equation for memory usage 'x' is the sum of the memory for each variable.
memory_used_bytes = mem_n_bytes + mem_m_bytes + mem_output_bytes

# The final output should be in the format x:o.
# The numbers that form the final result (the "equation") are calculated above.
# x = memory_used_bytes
# o = output_area
print(f"{memory_used_bytes}:{output_area}")
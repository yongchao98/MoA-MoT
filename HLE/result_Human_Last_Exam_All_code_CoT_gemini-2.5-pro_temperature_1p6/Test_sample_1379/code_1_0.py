# Define the dimensions of the rectangle
N = 1000000001
M = 1000000001

# Calculate the maximal coverable area.
# In integer arithmetic, this is equivalent to floor((N * M) / 2) * 2.
# Python's arbitrary-precision integers handle the large product automatically.
max_area = (N * M // 2) * 2

# Calculate the memory usage for an optimal C program's variables.
# 'n' and 'm' can be stored in 32-bit signed integers (4 bytes each).
mem_n_bytes = 4
mem_m_bytes = 4

# The product N*M is ~10^18, which requires a 64-bit signed integer
# ('long long' in C, 8 bytes) to avoid overflow.
mem_output_bytes = 8

# Total memory for the specified variables.
total_memory = mem_n_bytes + mem_m_bytes + mem_output_bytes

# Print the final result in the format x:o
# where x is memory in bytes and o is the maximal area.
print(f"{total_memory}:{max_area}")
# Plan:
# 1. Define the dimensions N and M.
# 2. Calculate the maximal covered area, which is the largest even number <= N * M.
#    This is equivalent to (N * M) // 2 * 2.
# 3. Calculate the memory usage in bytes for an optimal C program.
#    - N (1000000001) fits in an 'int' (4 bytes).
#    - M (1000000001) fits in an 'int' (4 bytes).
#    - The output area (~10^18) requires a 'long long' (8 bytes).
#    - Total memory = 4 + 4 + 8 = 16 bytes.
# 4. Print the memory and the calculated area in the format "x:o".

# Dimensions of the rectangle
N = 1000000001
M = 1000000001

# Calculate the maximal covered area using integer arithmetic
max_area = (N * M) // 2 * 2

# Memory usage for variables (n, m, output) in an optimal C program
# sizeof(int) for n, sizeof(int) for m, sizeof(long long) for output
memory_in_bytes = 4 + 4 + 8

# Print the final result in the specified format "x:o"
print(f"{memory_in_bytes}:{max_area}")
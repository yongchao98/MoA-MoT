import sys

# Define the dimensions of the rectangle as per the problem statement.
N = 1000000001
M = 1000000001

# Part 1: Calculate the maximal covered area 'o'.
# Each brick has an area of 2. The maximum covered area is the largest
# even number less than or equal to the total area (N * M).
# In integer arithmetic, this is calculated as (N * M) // 2 * 2.
total_area = N * M
max_covered_area = (total_area // 2) * 2

print("--- Calculating the Maximal Covered Area (o) ---")
print(f"The equation to find the maximal covered area 'o' is: ({N} * {M}) // 2 * 2")
print(f"Resulting area 'o': {max_covered_area}\n")


# Part 2: Calculate the memory usage 'x' for an equivalent C program.
# The values for N, M, and the result are large (~10^18),
# so they would be stored in 'long long' variables in C.
# A 'long long' typically occupies 8 bytes.
# We need variables for n, m, and the output.
sizeof_long_long = 8  # in bytes
num_variables = 3     # n, m, output
memory_usage = num_variables * sizeof_long_long

print("--- Calculating Memory Usage (x) ---")
print(f"In C, N, M, and the output would be stored as 'long long'.")
print(f"The size of 'long long' is typically {sizeof_long_long} bytes.")
print(f"The equation to find the memory usage 'x' is: {num_variables} variables * {sizeof_long_long} bytes/variable")
print(f"Resulting memory 'x' (in bytes): {memory_usage}\n")


# Final Answer: Output the result in the format x:o
print("--- Final Answer ---")
print("The result in the format x:o (memory_in_bytes:max_covered_area) is:")
print(f"{memory_usage}:{max_covered_area}")

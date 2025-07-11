# Define the dimensions of the rectangle
n = 1000000001
m = 1000000001

# --- Calculate the maximal coverable area ---
# The area of one brick is 2. The total area covered must be an even number.
# The maximum area is the largest even number less than or equal to the total area (n * m).
# This can be calculated using integer division.
total_area = n * m
max_covered_area = (total_area // 2) * 2

# --- Calculate the memory usage for an optimal C program ---
# For n and m (1,000,000,001), a 32-bit integer is sufficient (4 bytes).
# For the output (~10^18), a 64-bit integer is required (8 bytes).
# Total memory = memory for n + memory for m + memory for output
memory_for_n = 4  # bytes for a 32-bit int
memory_for_m = 4  # bytes for a 32-bit int
memory_for_output = 8 # bytes for a 64-bit long long
total_memory_bytes = memory_for_n + memory_for_m + memory_for_output

# Print the result in the specified "x:o" format
print(f"{total_memory_bytes}:{max_covered_area}")

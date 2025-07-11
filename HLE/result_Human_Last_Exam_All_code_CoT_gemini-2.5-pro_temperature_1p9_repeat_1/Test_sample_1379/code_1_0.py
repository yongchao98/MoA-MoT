# Plan:
# 1. Define the dimensions N and M.
n = 1000000001
m = 1000000001

# 2. Calculate the memory usage for variables in a C program.
# We'd use long long for n (8 bytes), long long for m (8 bytes),
# and unsigned long long for the result (8 bytes).
# Total memory = 8 + 8 + 8 = 24 bytes.
memory_usage = 24

# 3. Calculate the maximal coverable area.
# Since n and m are both odd, their product is odd.
# The maximal area is n * m - 1.
maximal_area = n * m - 1

# 4. Print the result in the format x:o
# The equation is: maximal_area = 1000000001 * 1000000001 - 1
# The numbers in the equation are n, m, and the result, which are all used in this script.
print(f"{memory_usage}:{maximal_area}")

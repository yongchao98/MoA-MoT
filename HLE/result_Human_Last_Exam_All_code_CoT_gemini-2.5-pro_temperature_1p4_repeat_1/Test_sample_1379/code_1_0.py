# Set the dimensions of the rectangle as per the problem statement.
n = 1000000001
m = 1000000001

# Calculate the maximal area that can be covered by 2x1 bricks.
# The total area of the rectangle is n * m. Since each brick has an area of 2,
# the maximum coverable area is the largest even number not exceeding the total area.
# This can be computed using integer division.
maximal_area = (n * m // 2) * 2

# The problem asks for the memory usage (x) of an equivalent C program.
# In C, to handle numbers of this magnitude (n * m is ~10^18), we would use
# the `long long` data type, which is typically 8 bytes.
# We need variables for n, m, and the output.
# So, the memory usage is 3 variables * 8 bytes/variable.
memory_used_bytes = 24

# The final answer format is x:o, where x is memory and o is the output.
# The instruction "output each number in the final equation" means we should
# print the two components of our final answer string.
print(f"{memory_used_bytes}:{maximal_area}")
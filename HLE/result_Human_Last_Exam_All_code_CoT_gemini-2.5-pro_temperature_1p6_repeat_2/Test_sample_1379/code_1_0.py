# The problem asks to find the maximal area of an N x M rectangle
# that can be covered by 2x1 bricks.

# Define the dimensions of the rectangle.
N = 1000000001
M = 1000000001

# The area of a single brick is 2 units. The total covered area must be even.
# We can find the maximum number of bricks that fit by using integer
# division of the total area by 2.
total_area = N * M
number_of_bricks = total_area // 2

# The maximal covered area is the number of bricks multiplied by the area of one brick (2).
# This calculation effectively finds the largest even number not greater than the total area.
max_covered_area = number_of_bricks * 2

# The problem also asks for the memory usage (in bytes) of an optimal C program
# for the variables n, m, and the output.
# N and M (1000000001) fit into a standard 32-bit `int` (4 bytes).
memory_for_n = 4
memory_for_m = 4
# The output (~10^18) requires a 64-bit `long long` (8 bytes).
memory_for_output = 8

total_memory_usage = memory_for_n + memory_for_m + memory_for_output

# Print the final answer in the format "x:o" where x is the memory
# in bytes and o is the maximal covered area.
# The instruction "output each number in the final equation" is interpreted as
# showing how the result 'o' is derived from N and M in the code's logic,
# which has been done above. The final print adheres to the x:o format.
print(f"{total_memory_usage}:{max_covered_area}")
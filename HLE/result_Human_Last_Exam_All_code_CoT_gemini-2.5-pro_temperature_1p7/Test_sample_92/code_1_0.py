# The bins where the marble escapes or melts.
escape_bin = 2025
melt_bin = 2024

# The starting position of the marble.
start_bin = 0

# The probability of escaping from bin n, p_n, follows the relation p_n = A*n + B.
# We have a system of two linear equations from the boundary conditions:
# p_2024 = 0  => A * 2024 + B = 0
# p_2025 = 1  => A * 2025 + B = 1

# We can solve this system for A and B.
# (A * 2025 + B) - (A * 2024 + B) = 1 - 0  => A = 1
A = 1
# Substitute A=1 into the first equation: 1 * 2024 + B = 0 => B = -2024
B = -melt_bin

# The equation for the escape probability from bin n is p_n = n - 2024.
# We want to find the probability of escaping from the start_bin.
prob_escape_from_start = A * start_bin + B

# We need to print each number in the final equation.
print("The probability of escaping from bin n is given by the equation: p_n = A*n + B")
print(f"From the boundary conditions p_{melt_bin} = 0 and p_{escape_bin} = 1, we solve for A and B.")
print(f"A = {A}")
print(f"B = {B}")
print(f"So, p_n = {A}*n + ({B}) = n - {-B}")
print(f"For the starting bin n = {start_bin}:")
print(f"p_{start_bin} = {A}*{start_bin} - {-B}")
print(f"p_{start_bin} = {prob_escape_from_start}")

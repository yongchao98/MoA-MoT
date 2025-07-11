# Define the positions of the escape portal and the torch
escape_bin = 2025
melt_bin = 2024

# The starting position of the marble
start_bin = 0

# For this type of unbiased random walk, the probability p_n of reaching
# 'escape_bin' before 'melt_bin' starting from 'n' is a linear function:
# p_n = A*n + B.

# We set up a system of linear equations using the boundary conditions:
# p_escape = 1  => A * escape_bin + B = 1
# p_melt = 0    => A * melt_bin + B = 0

# We can solve this system for A and B.
# (A * escape_bin + B) - (A * melt_bin + B) = 1 - 0
# A * (escape_bin - melt_bin) = 1
A = 1 / (escape_bin - melt_bin)

# Substitute A back into the second equation:
# A * melt_bin + B = 0
# B = -A * melt_bin
B = -A * melt_bin

# The solution is p_n = A*n + B.
# We need to find the probability for the starting bin, n = 0.
p_start = A * start_bin + B

# Let's print the full equation for clarity.
# p_n = A*n + B
# p_0 = A*0 + B = B
# B = -A * melt_bin = -(1 / (escape_bin - melt_bin)) * melt_bin
# B = -melt_bin / (escape_bin - melt_bin)
numerator = -melt_bin
denominator = escape_bin - melt_bin
result = numerator / denominator

print(f"The probability is given by the formula (n - melt_bin) / (escape_bin - melt_bin).")
print(f"For a start at n = {start_bin}:")
print(f"p_{start_bin} = ({start_bin} - {melt_bin}) / ({escape_bin} - {melt_bin})")
print(f"p_{start_bin} = ({start_bin - melt_bin}) / ({escape_bin - melt_bin})")
print(f"p_{start_bin} = {p_start}")

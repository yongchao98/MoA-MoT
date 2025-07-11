import collections

# Step 1: Calculate N_main, the number of involutions on 8 elements.
# The number of involutions I_n on a set of n elements follows the recurrence relation:
# I_n = I_{n-1} + (n-1) * I_{n-2}, with base cases I_0 = 1, I_1 = 1.
I = collections.defaultdict(int)
I[0] = 1
I[1] = 1
for n in range(2, 9):
    I[n] = I[n-1] + (n - 1) * I[n-2]

N_main = I[8]

# Step 2: Calculate N_anti.
# The number of configurations symmetric about the anti-diagonal is also I_8.
N_anti = N_main

# Step 3: Calculate N_both, the number of configurations symmetric about both diagonals.
# Let a(2m) be this number for a 2m x 2m board. The recurrence relation is:
# a(2m) = 2 * a(2m-2) + (2m-2) * a(2m-4), with base cases a(0) = 1, a(2) = 2.
# We need to find a(8), which corresponds to m=4.
# Let's define b[m] = a(2m).
b = collections.defaultdict(int)
b[0] = 1  # Corresponds to a(0) for a 0x0 board
b[1] = 2  # Corresponds to a(2) for a 2x2 board
for m in range(2, 5): # m=2,3,4
    b[m] = 2 * b[m - 1] + (2 * m - 2) * b[m - 2]

N_both = b[4]

# Step 4: Apply the Principle of Inclusion-Exclusion.
# Total = N_main + N_anti - N_both
total_configs = N_main + N_anti - N_both

# Print the final equation
print(f"The number of configurations symmetric along the main diagonal is {N_main}.")
print(f"The number of configurations symmetric along the anti-diagonal is {N_anti}.")
print(f"The number of configurations symmetric along both diagonals is {N_both}.")
print(f"Using the Principle of Inclusion-Exclusion, the total number of possible configurations is:")
print(f"{N_main} + {N_anti} - {N_both} = {total_configs}")

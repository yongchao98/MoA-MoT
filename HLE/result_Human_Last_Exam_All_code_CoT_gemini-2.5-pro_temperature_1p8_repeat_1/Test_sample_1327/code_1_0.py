# The given sequence
T = [2, 11, 23, 51, 119]

# 1. Calculate the sequence of differences between terms
d = []
for i in range(len(T) - 1):
    diff = T[i+1] - T[i]
    d.append(diff)

# 2. Find the recursive pattern in the differences
# The pattern appears to be d[i+1] = 2 * d[i] + C[i] starting from i=1.
# Let's find the sequence C.
C = []
# Calculate C for d[2] = 2 * d[1] + C[0]
# 28 = 2 * 12 + c1 => c1 = 4
c1 = d[2] - 2 * d[1]
C.append(c1)
# Calculate C for d[3] = 2 * d[2] + C[1]
# 68 = 2 * 28 + c2 => c2 = 12
c2 = d[3] - 2 * d[2]
C.append(c2)

# 3. We notice a pattern in C: C[i+1] = 3 * C[i]
# Predict the next value in the C sequence
next_c = C[-1] * 3

# 4. Predict the next value in the difference sequence d
next_d = 2 * d[-1] + next_c

# 5. Predict the next value in the original sequence T
next_T = T[-1] + next_d

# Print the final calculation as requested
print(f"The next term is found by adding the next difference to the last term.")
print(f"The next difference is calculated to be {next_d}.")
print(f"So, the final equation is:")
print(f"{T[-1]} + {next_d} = {next_T}")
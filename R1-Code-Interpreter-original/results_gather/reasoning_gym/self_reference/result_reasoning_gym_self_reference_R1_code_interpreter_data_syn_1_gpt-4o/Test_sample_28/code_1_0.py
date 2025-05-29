# Possible values for T and F
possible_solutions = []

# Check each possible combination of T and F
for T in range(8):
    F = 7 - T
    # Check if the combination satisfies all conditions
    if (T == 3 and F == 4):
        possible_solutions.append((T, F))

# Output the number of possible solutions
print(len(possible_solutions))
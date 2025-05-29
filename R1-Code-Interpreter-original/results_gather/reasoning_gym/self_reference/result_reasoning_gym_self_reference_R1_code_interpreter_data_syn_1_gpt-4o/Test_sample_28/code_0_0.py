# Possible values for T and F
possible_solutions = []

# Check each possible combination of T and F
for T in range(8):
    F = 7 - T
    # Check if the combination satisfies all conditions
    if (T in [5, 6, 7]) and (F in [0, 1, 2]) and (T == 3 or F == 1) and (T in [2, 3, 5, 7]) and (F in [4, 6]):
        possible_solutions.append((T, F))

# Output the number of possible solutions
print(len(possible_solutions))
# Possible values for T (true statements) and F (false statements)
possible_T = [2, 3, 5, 7]  # Prime numbers
possible_F = [4, 6]        # Composite numbers

# Count the number of valid solutions
valid_solutions = 0

# Check each combination of T and F
for T in possible_T:
    F = 7 - T
    if F in possible_F:
        # Check the conditions for Statement 5
        if (T == 3 and F != 2) or (F == 2 and T != 3):
            valid_solutions += 1

print(valid_solutions)
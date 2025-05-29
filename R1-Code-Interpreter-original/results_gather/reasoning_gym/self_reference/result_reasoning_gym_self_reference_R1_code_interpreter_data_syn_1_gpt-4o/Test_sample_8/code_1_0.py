# Check the number of possible solutions
possible_solutions = 0

# Iterate over possible values of T (number of true statements)
for T in range(2, 8):
    F = 7 - T  # Calculate the number of false statements
    # Check if the conditions are satisfied
    if (T >= 2 and T >= 3 and T != 7 and F != 3 and
        (T in [2, 3, 5, 7]) and (F in [4, 6])):
        possible_solutions += 1

print(possible_solutions)
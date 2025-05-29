# Possible values for T and F
possible_T = [4, 5, 6, 7]
possible_F = [0, 1, 2, 3, 4, 5]

# Prime numbers for T
prime_T = [2, 3, 5, 7]

# Composite numbers for F
composite_F = [4, 6, 8, 9, 10, 12, 14]

# Check possible solutions
solutions = 0
for T in possible_T:
    F = 7 - T
    if F in possible_F and T in prime_T and F in composite_F:
        # Check Statement 5 condition
        if (T == 4 and F != 3) or (F == 3 and T != 4):
            solutions += 1

print(solutions)
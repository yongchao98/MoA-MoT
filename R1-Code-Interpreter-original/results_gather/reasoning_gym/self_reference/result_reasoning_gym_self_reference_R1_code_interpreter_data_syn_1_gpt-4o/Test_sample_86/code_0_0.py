# Possible values for T and F
possible_solutions = 0

# Check each possible value of T
for T in range(1, 8):
    F = 7 - T
    # Check if T is prime
    if T in [2, 3, 5, 7]:
        # Check if F is composite
        if F in [4, 6, 8]:
            possible_solutions += 1

print(possible_solutions)
from itertools import permutations

# Given numbers
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Define the equations as functions
def check_constraints(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = values
    return (
        D == 3 * L and
        F + L == 17 and
        I + J == 33 and
        I - D == -17 and
        I - E == 21 and
        D - C == 29 and
        A == 4 * B and
        F - C == -14 and
        L == 3 * J and
        L - D == -30 and
        F - B == -22
    )

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    if check_constraints(perm):
        # If a valid permutation is found, print it
        print(f"<<<{list(perm)}>>>")
        break
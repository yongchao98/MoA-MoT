from itertools import permutations

# Given numbers
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Constraints as functions
def constraints(assignment):
    A, B, C, D, E, F, G, H, I, J, K, L, M = assignment
    return (
        F + J == 15 and
        F - G == -91 and
        D == 3.5 * K and
        E + G == 141 and
        A - K == 13 and
        F + M == 29 and
        D == 1.4 * F and
        L == 4.0 * D and
        C == 3.6 * J and
        A + M == 39 and
        C == 2.4 * A
    )

# Try all permutations of the numbers
for perm in permutations(numbers):
    if constraints(perm):
        # If a valid assignment is found, print it
        print(f"<<<{list(perm)}>>>")
        break
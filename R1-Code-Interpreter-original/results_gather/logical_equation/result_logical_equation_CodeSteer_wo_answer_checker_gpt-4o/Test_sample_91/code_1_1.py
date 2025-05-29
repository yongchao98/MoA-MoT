from itertools import permutations

# Given numbers
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Constraints as functions
def check_constraints(assignment):
    A, B, C, D, E, F, G, H, I, J, K, L, M = assignment
    return (
        I == 2.4 * F and
        D == 4.8 * C and
        D == 1.6 * F and
        B == 1.6 * L and
        M == 3.5 * E and
        G - J == -17 and
        H - B == 34 and
        B == 3.2 * C and
        H > E and
        H + L == 60 and
        I - C == 31 and
        H > A
    )

# Iterate over all permutations and check constraints
for perm in permutations(numbers):
    if check_constraints(perm):
        print(f"<<<{list(perm)}>>>")
        break
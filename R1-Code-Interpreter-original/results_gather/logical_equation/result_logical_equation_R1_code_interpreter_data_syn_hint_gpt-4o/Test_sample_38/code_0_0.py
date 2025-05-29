from itertools import permutations

# Given numbers
numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Constraints
def check_constraints(values):
    A, B, C, D, E, F, G, H, I, J, K = values
    return (
        C + H == 124 and
        B + F == 305 and
        G == 3 * J and
        H - G == -45 and
        C + J == 137 and
        J - A == -11 and
        E - J == 59 and
        F == 3 * E and
        A + H == 30 and
        A - C == -94
    )

# Find the correct permutation
for perm in permutations(numbers):
    if check_constraints(perm):
        print(list(perm))
        break
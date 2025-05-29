from itertools import permutations

# Given numbers
numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Constraints
def satisfies_constraints(assignment):
    A, B, C, D, E, F, G, H, I, J, K = assignment
    return (
        K > A and
        H == 2.4 * A and
        H == 3.0 * F and
        B - I == -216 and
        C - J == -72 and
        I == 1.5 * G and
        I == 3.0 * J and
        B - H == -39 and
        D - I == -104 and
        J - F == 59 and
        B - F == -7
    )

# Try all permutations
for perm in permutations(numbers):
    if satisfies_constraints(perm):
        print(f"<<<{list(perm)}>>>")
        break
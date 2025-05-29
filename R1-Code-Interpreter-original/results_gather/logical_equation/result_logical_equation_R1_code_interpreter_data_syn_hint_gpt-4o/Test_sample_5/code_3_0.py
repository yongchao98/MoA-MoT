from itertools import permutations

# Given numbers
numbers = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Constraints
def satisfies_constraints(perm):
    A, B, C, D, E, F, G, H, I = perm
    return (
        G == 4 * A and
        G > A and
        C == 4 * D and
        B + E == 148 and
        A - E == -117 and
        D - F == -71 and
        C + H == 39
    )

# Find the correct permutation
for perm in permutations(numbers):
    if satisfies_constraints(perm):
        result = list(perm)
        break

print(result)
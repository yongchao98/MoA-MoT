from itertools import permutations

# Given numbers
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Constraints as functions
def check_constraints(perm):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    return (
        F + J == 86 and
        B == 4.5 * A and
        A + K == 15 and
        M == 3.0 * L and
        H == 1.5 * A and
        J > D and
        A - M == 7 and
        B + E == 52 and
        A > G and
        H > L
    )

# Iterate over all permutations
for perm in permutations(numbers):
    if check_constraints(perm):
        # Output the result in alphabetical order
        result = list(perm)
        print(f"<<<{result}>>>")
        break
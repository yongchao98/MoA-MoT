from itertools import permutations

# Given numbers
numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Constraints
def check_constraints(values):
    A, B, C, D, E, F, G, H, I, J, K = values
    return (
        H + J == 130 and
        C + I == 177 and
        C + G == 75 and
        G == 3 * F and
        B - E == 150 and
        J - C == -18 and
        E + G == 123 and
        D + F == 96 and
        E + K == 78 and
        H + K == 124
    )

# Find the correct permutation
for perm in permutations(numbers):
    if check_constraints(perm):
        print(list(perm))
        break
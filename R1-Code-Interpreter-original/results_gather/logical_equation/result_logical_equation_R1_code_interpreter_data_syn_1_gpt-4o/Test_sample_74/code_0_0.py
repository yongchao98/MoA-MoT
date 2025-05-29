from itertools import permutations

# Possible values for the letters
values = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Constraints
def check_constraints(assignment):
    A, B, C, D, E, F, G, H, I, J, K, L, M = assignment
    return (
        E + F == 8 and
        C - D == 14 and
        A + M == 35 and
        D + K == 3 and
        G + J == 60 and
        B > F and
        H + I == 60 and
        A - K == 5 and
        G - H == 26 and
        G - F == 45
    )

# Find the correct permutation
for perm in permutations(values):
    if check_constraints(perm):
        print(list(perm))
        break
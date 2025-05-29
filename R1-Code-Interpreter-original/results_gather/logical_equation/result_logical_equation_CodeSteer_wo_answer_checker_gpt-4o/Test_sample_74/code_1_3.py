from itertools import permutations

# Define the possible values
possible_values = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Define the constraints as functions
def satisfies_constraints(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = values
    return (
        E + F == 8 and
        C - D == 14 and
        A + M == 35 and
        D + K == 3 and
        G + J == 60 and
        H + I == 60 and
        A - K == 5 and
        G - H == 26 and
        G - F == 45 and
        B > F
    )

# Iterate over all permutations of the possible values
for perm in permutations(possible_values):
    if satisfies_constraints(perm):
        result = list(perm)
        print(f"<<<{result}>>>")
        break
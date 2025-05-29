from itertools import permutations

# Possible values for the letters
values = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Define the constraints
def satisfies_constraints(assignment):
    A, B, C, D, E, F, G, H, I, J, K = assignment
    return (
        D + J == 252 and
        D == 3 * K and
        H - A == 44 and
        J > E and
        C == 2 * K and
        D > A and
        A + F == 52 and
        D + E == 228 and
        B == 3 * F and
        E + I == 124
    )

# Iterate over all permutations of the values
for perm in permutations(values):
    if satisfies_constraints(perm):
        print(f"<<<{list(perm)}>>>")
        break
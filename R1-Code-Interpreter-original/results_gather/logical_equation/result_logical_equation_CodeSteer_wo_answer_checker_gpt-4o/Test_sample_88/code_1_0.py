from itertools import permutations

# List of possible values
values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Function to check if a permutation satisfies all constraints
def satisfies_constraints(perm):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    return (
        J == 3.5 * D and
        K - F == -12 and
        G + M == 26 and
        J + K == 31 and
        C - D == 13 and
        G == 1.6 * M and
        L > G and
        A + L == 124 and
        M - A == -18 and
        C == 3.0 * H and
        J > D
    )

# Iterate over all permutations of the values
for perm in permutations(values):
    if satisfies_constraints(perm):
        # If all constraints are satisfied, print the solution
        print(f"<<<{list(perm)}>>>")
        break
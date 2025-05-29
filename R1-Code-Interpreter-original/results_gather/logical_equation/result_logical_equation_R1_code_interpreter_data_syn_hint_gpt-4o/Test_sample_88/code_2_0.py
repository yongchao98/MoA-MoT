from itertools import permutations

# Possible values for the letters
values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all equations
    if (J == 3.5 * D and
        K - F == -12 and
        G + M == 26 and
        J + K == 31 and
        C - D == 13 and
        G == 1.6 * M and
        A + L == 124 and
        M - A == -18 and
        C == 3.0 * H):
        
        # Check inequalities
        if J > D and L > G:
            # If all constraints are satisfied, print the solution
            print([A, B, C, D, E, F, G, H, I, J, K, L, M])
            break
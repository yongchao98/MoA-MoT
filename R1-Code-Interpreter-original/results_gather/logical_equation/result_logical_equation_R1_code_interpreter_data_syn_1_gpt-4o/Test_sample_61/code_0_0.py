from itertools import permutations

# Possible values for the letters
values = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all the conditions
    if (C + M == 43 and
        C == 1.5 * L and
        L == 2.0 * J and
        A + K == 81 and
        C - B == -35 and
        H - K == -20 and
        F == 4.8 * J and
        K - M == 8 and
        D + E == 4 and
        K == 2.4 * C):
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
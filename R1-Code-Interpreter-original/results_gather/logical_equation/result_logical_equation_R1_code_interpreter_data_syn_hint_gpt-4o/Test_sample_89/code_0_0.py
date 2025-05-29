from itertools import permutations

# Possible values for the letters
values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all conditions
    if (A - C == 26 and
        M > L and
        C == 2.0 * D and
        E + I == 57 and
        L == 1.6 * C and
        G - K == -12 and
        K == 1.5 * C and
        A == 3.6 * C and
        M == 4.8 * D and
        C + F == 106 and
        H - G == 42):
        
        # If all conditions are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
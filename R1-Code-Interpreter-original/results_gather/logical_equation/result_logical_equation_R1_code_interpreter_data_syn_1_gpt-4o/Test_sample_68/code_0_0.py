from itertools import permutations

# Possible values for the letters
values = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all conditions
    if (B + L == 35 and
        J > D and
        L - J == -8 and
        I == 3.6 * H and
        B + I == 64 and
        H + I == 46 and
        C + L == 31 and
        G + K == 19 and
        M - J == 30 and
        K == 3.2 * E):
        
        # If all conditions are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
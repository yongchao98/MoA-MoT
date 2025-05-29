from itertools import permutations

# Possible values for the letters
values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all conditions
    if (L - I == -23 and
        I + M == 44 and
        D + K == 98 and
        F + I == 73 and
        F > J and
        B - C == 26 and
        B == 2.4 * A and
        C + F == 55 and
        F - J == 38 and
        H - F == 5 and
        J == 1.4 * L):
        
        # If all conditions are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
from itertools import permutations

# Possible values for the letters
values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all constraints
    if (L == 1.4 * K and
        B + C == 73 and
        J - G == -86 and
        D + H == 53 and
        B + K == 50 and
        L - A == -29 and
        F - I == -13 and
        J == 2.0 * K and
        C + M == 52 and
        I > L and
        D + M == 74 and
        M - B == -21):
        
        # If all constraints are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
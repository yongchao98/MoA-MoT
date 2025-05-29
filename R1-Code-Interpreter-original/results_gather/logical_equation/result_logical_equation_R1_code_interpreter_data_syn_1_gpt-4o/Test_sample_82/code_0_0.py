from itertools import permutations

# Possible values for the letters
values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all conditions
    if (C - E == -21 and
        J == 4.0 * C and
        H + K == 23 and
        D + L == 17 and
        C - A == -26 and
        G + M == 8 and
        F + G == 13 and
        H - G == 4 and
        A - C == 26 and
        I - J == -68 and
        F + M == 15):
        
        # If all conditions are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
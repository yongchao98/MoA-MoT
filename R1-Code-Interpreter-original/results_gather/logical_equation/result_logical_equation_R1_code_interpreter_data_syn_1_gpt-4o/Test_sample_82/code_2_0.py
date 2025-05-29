from itertools import permutations

# Possible values for the letters
values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over all possible values for C
for C in values:
    J = 4 * C
    A = C + 26
    I = J - 68
    
    # Check if J, A, and I are in the list of values
    if J in values and A in values and I in values:
        # Remove C, J, A, and I from the list of available values
        remaining_values = [v for v in values if v not in {C, J, A, I}]
        
        # Iterate over permutations of the remaining values
        for perm in permutations(remaining_values):
            D, E, F, G, H, K, L, M = perm
            
            # Check all conditions
            if (C - E == -21 and
                H + K == 23 and
                D + L == 17 and
                G + M == 8 and
                F + G == 13 and
                H - G == 4 and
                F + M == 15):
                
                # If all conditions are satisfied, print the result
                print([A, B, C, D, E, F, G, H, I, J, K, L, M])
                break
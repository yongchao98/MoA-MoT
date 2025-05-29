from itertools import permutations

# Possible values for the letters
values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Precompute valid (F, J) pairs
valid_FJ = [(F, J) for F in values for J in values if F - J == 38]

# Precompute valid (F, H) pairs
valid_FH = [(F, H) for F in values for H in values if H - F == 5]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check if (F, J) and (F, H) are valid pairs
    if (F, J) in valid_FJ and (F, H) in valid_FH:
        # Calculate dependent variables
        L = I - 23
        M = 44 - I
        C = 55 - F
        
        # Check all conditions
        if (L in values and M in values and C in values and
            D + K == 98 and
            F > J and
            B - C == 26 and
            B == 2.4 * A and
            J == 1.4 * L):
            
            # If all conditions are satisfied, print the result
            print([A, B, C, D, E, F, G, H, I, J, K, L, M])
            break
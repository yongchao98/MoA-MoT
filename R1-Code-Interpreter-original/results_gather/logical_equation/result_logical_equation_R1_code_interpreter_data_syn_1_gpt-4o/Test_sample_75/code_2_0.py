from itertools import permutations

# Possible values for the letters
values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Precompute valid (A, B) pairs
valid_AB = [(A, B) for A in values for B in values if B == 2.4 * A]

# Precompute valid (L, J) pairs
valid_LJ = [(L, J) for L in values for J in values if J == 1.4 * L]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check if (A, B) and (L, J) are valid pairs
    if (A, B) in valid_AB and (L, J) in valid_LJ:
        # Calculate dependent variables
        L = I - 23
        M = 44 - I
        F = 73 - I
        C = 55 - F
        J = F - 38
        H = F + 5
        
        # Check all conditions
        if (L in values and M in values and C in values and J in values and H in values and
            D + K == 98 and
            F > J and
            B - C == 26):
            
            # If all conditions are satisfied, print the result
            print([A, B, C, D, E, F, G, H, I, J, K, L, M])
            break
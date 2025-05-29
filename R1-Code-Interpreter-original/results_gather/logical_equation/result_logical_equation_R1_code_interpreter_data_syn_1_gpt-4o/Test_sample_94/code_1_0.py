from itertools import permutations

# Possible values for the letters
values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Pre-calculate possible pairs for G and I
possible_g_i_pairs = [(g, i) for g in values for i in values if g + i == 55]

# Iterate over possible G and I pairs
for G, I in possible_g_i_pairs:
    E = G - 5
    L = 57 - G
    F = I + 91
    
    # Check if E, L, F are in the list of values
    if E in values and L in values and F in values:
        # Remove G, I, E, L, F from the list of values
        remaining_values = [v for v in values if v not in {G, I, E, L, F}]
        
        # Iterate over permutations of the remaining values
        for perm in permutations(remaining_values):
            A, B, C, D, H, J, K, M = perm
            
            # Check all the conditions
            if (E > D and
                C + H == 27 and
                A - I == 5 and
                B - C == -9 and
                K > C and
                H + L == 10 and
                J - D == 14 and
                B - F == -81):
                
                # If all conditions are satisfied, print the result
                print([A, B, C, D, E, F, G, H, I, J, K, L, M])
                break
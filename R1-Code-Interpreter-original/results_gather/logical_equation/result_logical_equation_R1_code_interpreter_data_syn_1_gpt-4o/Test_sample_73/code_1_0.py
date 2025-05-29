from itertools import permutations

# Possible values for the letters
values = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Pre-calculate possible values for some variables based on constraints
possible_L = [v for v in values if 3 * v in values]  # L such that M = 3L
possible_H = [v for v in values if 2.5 * v in values]  # H such that F = 2.5H

# Iterate over permutations of the remaining values
for L in possible_L:
    M = 3 * L
    if M not in values:
        continue
    
    for H in possible_H:
        F = 2.5 * H
        if F not in values:
            continue
        
        for perm in permutations(values):
            A, B, C, D, E, G, I, J, K = perm
            
            # Calculate dependent variables
            if E - L == 9 and A + I == 60 and E + J == 17 and G - I == 21 and G - A == -9 and E + G == 46 and F + M == 8 and H + I == 17:
                # If all constraints are satisfied, print the result
                print([A, B, C, D, E, F, G, H, I, J, K, L, M])
                break
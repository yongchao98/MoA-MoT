from itertools import permutations

# Possible values for the letters
values = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Pre-calculate possible values for I and J based on I = 4.0J
possible_IJ = [(i, j) for i in values for j in values if i == 4 * j]

# Pre-calculate possible values for L, A, and F based on L = 3.2A and L = 1.6F
possible_LAF = [(l, a, f) for l in values for a in values for f in values if l == 3.2 * a and l == 1.6 * f]

# Pre-calculate possible values for D and H based on D = 1.5H
possible_DH = [(d, h) for d in values for h in values if d == 1.5 * h]

# Iterate over all pre-calculated possibilities
for (I, J) in possible_IJ:
    for (L, A, F) in possible_LAF:
        for (D, H) in possible_DH:
            # Calculate E based on A - E = 4
            E = A - 4
            if E not in values:
                continue
            
            # Calculate M based on H + M = 26
            M = 26 - H
            if M not in values:
                continue
            
            # Check remaining conditions
            if (D - B == -14 and
                F + L == 26 and
                H - I == -4):
                
                # Find B such that D - B = -14
                B = D + 14
                if B not in values:
                    continue
                
                # Ensure all values are unique
                if len(set([A, B, C, D, E, F, G, H, I, J, K, L, M])) == 13:
                    print([A, B, C, D, E, F, G, H, I, J, K, L, M])
                    break
from itertools import permutations

# Possible values for the letters
values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Apply constraints to reduce search space
    if I == 1.5 * D and F == 4.0 * I and M == 2.0 * E and K - D == -13 and D - H == -34:
        # Check all conditions
        if (B + I == 31 and
            K - E == -2 and
            B + H == 57 and
            C + F == 132 and
            A + G == 43 and
            M > B and
            I + L == 69):
            print([A, B, C, D, E, F, G, H, I, J, K, L, M])
            break
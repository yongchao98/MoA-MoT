from itertools import permutations

# Possible values for the letters
values = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all conditions
    if (C - D == -94 and
        C - F == 18 and
        B > E and
        A > E and
        G == 4 * E and
        B > F and
        C == 3 * F):
        print([A, B, C, D, E, F, G, H, I])
        break
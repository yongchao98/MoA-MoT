from itertools import permutations

# Possible values for the letters
values = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all the conditions
    if (B + I == 156 and
        B + D == 123 and
        D - A == -29 and
        G == 2.0 * D and
        B == 1.5 * C and
        A > D and
        A - I == -4 and
        C - F == 71):
        print([A, B, C, D, E, F, G, H, I])
        break
from itertools import permutations

# Possible values for the letters
values = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all constraints
    if (D - E == 24 and
        A - F == 35 and
        D == 3 * B and
        I == 4 * F and
        B == 3 * E and
        E + H == 124 and
        C > E):
        print([A, B, C, D, E, F, G, H, I])
        break
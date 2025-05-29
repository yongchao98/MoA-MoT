from itertools import permutations

# Possible values for the letters
values = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all conditions
    if (A + E == 107 and
        D + H == 137 and
        H > C and
        A + B == 63 and
        E - F == 79 and
        B + E == 116 and
        B - D == -85):
        print([A, B, C, D, E, F, G, H, I])
        break
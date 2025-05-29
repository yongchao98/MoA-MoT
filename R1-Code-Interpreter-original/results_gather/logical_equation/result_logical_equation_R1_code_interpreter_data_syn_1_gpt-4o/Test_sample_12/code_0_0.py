from itertools import permutations

# Possible values for the letters
values = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all constraints
    if (D + H == 36 and
        G - I == -79 and
        C - B == -20 and
        D == 3 * E and
        B == 4 * D and
        I > C and
        B + C == 52):
        print([A, B, C, D, E, F, G, H, I])
        break
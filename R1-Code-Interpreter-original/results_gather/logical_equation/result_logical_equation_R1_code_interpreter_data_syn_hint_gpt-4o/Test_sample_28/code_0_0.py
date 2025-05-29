from itertools import permutations

# Possible values for the letters
values = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all constraints
    if (G == 1.5 * C and
        E == 1.8 * D and
        D + I == 141 and
        B - H == -88 and
        H + I == 241 and
        H - G == 111 and
        H - B == 88 and
        A + I == 201):
        print([A, B, C, D, E, F, G, H, I])
        break
from itertools import permutations

# Possible values for the letters
values = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all conditions
    if (A + F == 141 and
        A + G == 124 and
        I - C == 26 and
        H - G == 33 and
        C == 2 * G and
        D - H == -27 and
        D - G == 6 and
        H - E == -84):
        print([A, B, C, D, E, F, G, H, I])
        break
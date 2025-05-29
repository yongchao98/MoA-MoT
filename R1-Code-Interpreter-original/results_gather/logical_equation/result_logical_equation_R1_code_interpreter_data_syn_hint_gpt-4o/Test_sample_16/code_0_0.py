from itertools import permutations

# Possible values for the letters
values = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all conditions
    if (C == 3.0 * D and
        C == 1.5 * F and
        D - C == -6 and
        F == 2.0 * D and
        D + E == 124 and
        A + H == 112 and
        I == 1.5 * A and
        D + G == 23):
        print([A, B, C, D, E, F, G, H, I])
        break
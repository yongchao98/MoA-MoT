from itertools import permutations

# Possible values for the letters
values = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all the conditions
    if (D == 3 * B and
        A + F == 52 and
        D - H == -71 and
        B + I == 30 and
        D + E == 13 and
        I == 3 * D and
        F == 4 * E):
        print([A, B, C, D, E, F, G, H, I])
        break
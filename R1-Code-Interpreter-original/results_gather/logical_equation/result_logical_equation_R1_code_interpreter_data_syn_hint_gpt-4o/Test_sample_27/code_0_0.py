from itertools import permutations

# Possible values for the letters
values = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all constraints
    if (C == 3 * H and
        A - I == 16 and
        I - H == 17 and
        B == 4 * I and
        E + F == 127 and
        F - G == -114 and
        A == 4 * C and
        B == 2.5 * D):
        print([A, B, C, D, E, F, G, H, I])
        break
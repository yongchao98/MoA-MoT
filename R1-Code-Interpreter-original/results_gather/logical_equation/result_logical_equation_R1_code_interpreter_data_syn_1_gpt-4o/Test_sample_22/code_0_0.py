from itertools import permutations

# Possible values for the letters
values = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all conditions
    if (B == 2.0 * G and
        I == 2.5 * F and
        B + I == 86 and
        A + G == 12 and
        I > F and
        F + I == 112 and
        C + H == 140 and
        A - H == -11):
        print([A, B, C, D, E, F, G, H, I])
        break
from itertools import permutations

# Possible values for the letters
values = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all conditions
    if (A - C == -44 and
        C > H and
        I == 3.0 * D and
        B - C == -79 and
        D == 3.0 * B and
        A == 4.0 * I and
        B - E == -3):
        print([A, B, C, D, E, F, G, H, I])
        break
from itertools import permutations

# Possible values for the letters
values = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all conditions
    if (G - E == -3 and
        D - C == -41 and
        C > B and
        F == 4.0 * H and
        I - G == 2 and
        A + E == 31 and
        A - B == 11):
        print([A, B, C, D, E, F, G, H, I])
        break
from itertools import permutations

# Possible values for the letters
values = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all the inequalities
    if (A + F == 196 and
        G - D == -6 and
        G - A == -72 and
        E == 3.0 * A and
        B > I and
        K - A == 75 and
        H - F == -41 and
        B - G == 36 and
        A - H == -5 and
        E - A == 150 and
        H - J == 64):
        print([A, B, C, D, E, F, G, H, I, J, K])
        break
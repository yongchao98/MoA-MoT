from itertools import permutations

# Possible values for the letters
values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all conditions
    if (B + I == 31 and
        I == 1.5 * D and
        K - E == -2 and
        B + H == 57 and
        C + F == 132 and
        K - D == -13 and
        M == 2.0 * E and
        A + G == 43 and
        M > B and
        D - H == -34 and
        F == 4.0 * I and
        I + L == 69):
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
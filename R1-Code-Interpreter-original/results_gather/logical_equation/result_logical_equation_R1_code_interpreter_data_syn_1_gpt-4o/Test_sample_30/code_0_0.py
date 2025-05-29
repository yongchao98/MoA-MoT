from itertools import permutations

# Possible values for the letters
values = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all the conditions
    if (C + J == 186 and
        H == 3 * F and
        C - B == -44 and
        E + H == 234 and
        E + I == 130 and
        I - J == -29 and
        D + E == 57 and
        K - B == -64 and
        C + D == 84 and
        B - G == 53):
        print([A, B, C, D, E, F, G, H, I, J, K])
        break
from itertools import permutations

# Possible values for the letters
values = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Iterate over all permutations of the values
for perm in permutations(values):
    B, C, A, D, E, F, G, H, I = perm
    
    # Check the constraints
    if C == 3 * B and A == 4 * C and F > A and D == 4 * H and E > B and A - H == 32 and F + H == 84:
        print([A, B, C, D, E, F, G, H, I])
        break
from itertools import permutations

# Possible values for the letters
remaining_values = [1, 3, 4, 16, 80, 121]

# Iterate over all permutations of the remaining values
for perm in permutations(remaining_values):
    A, B, C, D, F, G = perm
    E = 9
    H = 36
    I = 27
    
    # Check all constraints
    if (D == A + 77 and
        H > A and
        I > B and
        G == A + 2):
        print([A, B, C, D, E, F, G, H, I])
        break
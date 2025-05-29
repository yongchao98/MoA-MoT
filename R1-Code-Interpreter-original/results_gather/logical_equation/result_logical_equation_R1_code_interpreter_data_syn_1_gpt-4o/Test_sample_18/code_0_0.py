from itertools import permutations

# Possible values for the letters
values = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Iterate through all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all constraints
    if (A + C == 127 and
        D - E == -77 and
        A - B == 89 and
        B == 1.6 * F and
        E > I and
        E == 4.0 * F and
        I == 4.0 * H and
        H == 1.5 * C):
        
        # If all constraints are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I])
        break
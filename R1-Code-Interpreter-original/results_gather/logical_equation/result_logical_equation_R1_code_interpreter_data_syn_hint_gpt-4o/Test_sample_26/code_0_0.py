from itertools import permutations

# Possible values for the letters
values = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all the constraints
    if (B - I == -74 and
        C + E == 140 and
        A + G == 39 and
        E + H == 29 and
        B == 2 * A and
        E + I == 100 and
        A + F == 124 and
        A - H == -6):
        
        # If all constraints are satisfied, print the solution
        print([A, B, C, D, E, F, G, H, I])
        break
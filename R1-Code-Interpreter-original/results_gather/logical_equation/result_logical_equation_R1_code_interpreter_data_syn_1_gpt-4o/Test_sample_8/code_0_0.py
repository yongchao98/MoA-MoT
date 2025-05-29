from itertools import permutations

# Possible values for the letters
values = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check the constraints
    if B > F and D + E == 107 and D > F and D - F == 77 and C > B and C + I == 52 and B == 3.0 * F:
        # If all constraints are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I])
        break
from itertools import permutations

# Possible values for the letters
values = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all the equations
    if (C + I == 241 and
        E - H == -44 and
        I == 1.5 * H and
        D + E == 42 and
        D + H == 86 and
        F + G == 35 and
        D == 2.0 * F and
        A + I == 140):
        
        # If all equations are satisfied, print the result
        result = [A, B, C, D, E, F, G, H, I]
        print(result)
        break
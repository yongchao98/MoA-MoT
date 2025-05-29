from itertools import permutations

# Possible values for the letters
values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all the conditions
    if (G + I == 55 and
        E > D and
        G + L == 57 and
        G - E == 5 and
        C + H == 27 and
        A - I == 5 and
        F - I == 91 and
        B - C == -9 and
        K > C and
        H + L == 10 and
        J - D == 14 and
        B - F == -81):
        
        # If all conditions are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
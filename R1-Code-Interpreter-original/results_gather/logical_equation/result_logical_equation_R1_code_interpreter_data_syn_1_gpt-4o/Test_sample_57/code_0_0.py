from itertools import permutations

# Possible values for the letters
values = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all the conditions
    if (C - G == -36 and
        C - F == -17 and
        B + F == 68 and
        E + H == 137 and
        A + K == 375 and
        C + H == 124 and
        A > F and
        J > C and
        F - G == -19 and
        A > C and
        C + G == 42):
        
        # If all conditions are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I, J, K])
        break
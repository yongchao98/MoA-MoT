from itertools import permutations

# Possible values for the letters
values = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all the constraints
    if (M == 3 * L and
        F == 2.5 * H and
        E - L == 9 and
        A + I == 60 and
        E + J == 17 and
        G - I == 21 and
        G - A == -9 and
        E + G == 46 and
        F + M == 8 and
        H + I == 17):
        
        # If all constraints are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
from itertools import permutations

# Possible values for the letters
values = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all constraints
    if (F == 2 * E and
        B > E and
        A - H == 5 and
        D + L == 48 and
        B > J and
        B + J == 16 and
        A + B == 22 and
        A == 1.4 * E and
        G + J == 25 and
        C + E == 21):
        
        # If all constraints are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
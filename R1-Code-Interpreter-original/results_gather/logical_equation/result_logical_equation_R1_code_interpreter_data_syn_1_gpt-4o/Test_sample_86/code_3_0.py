from itertools import permutations

# Possible values for the letters
values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all constraints
    if (C + H == 52 and
        G + I == 13 and
        H == 2.4 * E and
        F + J == 57 and
        A - I == 25 and
        J - I == 4 and
        D + G == 34 and
        F + L == 95 and
        A == 4.0 * J and
        F + G == 60 and
        J == 1.4 * K):
        
        # If all constraints are satisfied, print the solution
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
from itertools import permutations

# Possible values for the letters
values = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all constraints
    if (A - G == -13 and
        B + C == 15 and
        H == 4.8 * B and
        I - B == -3 and
        L - J == -5 and
        B > A and
        I == 2.0 * F and
        C - H == -14 and
        B - I == 3 and
        A - F == 2):
        
        # If all constraints are satisfied, print the solution
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
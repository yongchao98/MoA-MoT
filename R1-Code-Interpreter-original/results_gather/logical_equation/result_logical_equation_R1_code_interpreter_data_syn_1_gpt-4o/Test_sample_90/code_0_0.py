from itertools import permutations

# Possible values for the letters
values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all constraints
    if (E + K == 81 and
        K > J and
        E == 1.5 * G and
        B + J == 5 and
        M - F == -6 and
        K == 4.5 * M and
        M > B and
        F == 3.2 * C and
        C - E == -31 and
        L == 1.4 * C and
        A + M == 38 and
        F + I == 112):
        
        # If all constraints are satisfied, print the solution
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
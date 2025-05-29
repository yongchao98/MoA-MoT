from itertools import permutations

# Possible values for the letters
values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all conditions
    if (D - A == 43 and
        F == 1.5 * C and
        A + F == 43 and
        C == 1.6 * L and
        L == 3.0 * J and
        F - J == 31 and
        H > L and
        M == 4.5 * G and
        E + G == 12 and
        J + K == 33 and
        F == 3.6 * G and
        M > A):
        
        # If all conditions are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
from itertools import permutations

# List of numbers to assign
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all the conditions
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
        
        # If all conditions are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
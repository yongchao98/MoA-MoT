from itertools import permutations

# Given numbers
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Iterate over all permutations of the numbers for the involved variables
for perm in permutations(numbers, 10):
    A, B, C, D, F, G, H, K, L, M = perm
    
    # Check all equations
    if (B == F - 21 and
        L == 73 - F and
        A == F + 22 and
        K == L - 42 and
        C == 43 - B and
        H == 60 - L and
        M == 17 - H and
        C == 2.4 * H and
        D == 4 - K and
        M == G / 2.5):
        
        # If all equations are satisfied, print the result
        result = [A, B, C, D, E, F, G, H, I, J, K, L, M]
        print(result)
        break
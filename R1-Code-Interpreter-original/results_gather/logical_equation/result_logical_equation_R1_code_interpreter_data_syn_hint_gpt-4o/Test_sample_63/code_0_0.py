from itertools import permutations

# Given numbers
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all conditions
    if (A - E == 4 and
        E - A == -4 and
        I == 4.0 * J and
        D - B == -14 and
        H + M == 26 and
        L == 3.2 * A and
        L == 1.6 * F and
        F + L == 26 and
        H - I == -4 and
        D == 1.5 * H):
        
        # If all conditions are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
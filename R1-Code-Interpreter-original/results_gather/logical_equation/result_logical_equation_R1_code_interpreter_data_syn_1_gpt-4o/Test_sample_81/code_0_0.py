from itertools import permutations

# Given numbers
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all conditions
    if (F + G == 46 and
        J - M == 80 and
        J - F == 86 and
        H - L == 26 and
        B + M == 31 and
        B == 1.5 * F and
        C == 4.0 * A and
        L == 1.5 * M and
        G == 3.6 * F and
        D == 1.5 * E and
        L - I == -21):
        
        # If all conditions are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
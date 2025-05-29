from itertools import permutations

# Given numbers
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all conditions
    if (D == 3 * I and
        C + H == 20 and
        A > F and
        F + G == 38 and
        L - G == -26 and
        C + F == 17 and
        E + L == 60 and
        M == 4.5 * L and
        A - K == -17 and
        D + H == 8):
        
        # If all conditions are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
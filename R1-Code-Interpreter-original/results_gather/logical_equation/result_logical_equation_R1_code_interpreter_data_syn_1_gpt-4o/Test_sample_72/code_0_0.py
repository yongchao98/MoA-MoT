from itertools import permutations

# Given numbers
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all conditions
    if (F + J == 86 and
        B == 4.5 * A and
        A + K == 15 and
        M == 3.0 * L and
        H == 1.5 * A and
        J > D and
        A - M == 7 and
        B + E == 52 and
        A > G and
        H > L):
        # If all conditions are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
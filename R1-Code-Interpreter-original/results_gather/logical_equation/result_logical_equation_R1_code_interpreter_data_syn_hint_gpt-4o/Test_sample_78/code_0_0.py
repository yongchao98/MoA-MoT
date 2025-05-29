from itertools import permutations

# Given numbers
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all conditions
    if (H == 2.4 * D and
        G - A == -9 and
        E == 4.0 * G and
        F - G == 8 and
        B + E == 78 and
        M - K == 40 and
        D + M == 55 and
        A == 3.2 * K and
        I - C == 94 and
        G + L == 43 and
        L - A == 20):
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
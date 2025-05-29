from itertools import permutations

# Given numbers
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all conditions
    if (A == 3.0 * L and
        A + M == 39 and
        L + M == 29 and
        B == 3.2 * L and
        A == 1.5 * I and
        C - B == 80 and
        E + K == 43 and
        C > D and
        E == 2.4 * A and
        C + D == 124 and
        H - A == -12 and
        F + G == 47):
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
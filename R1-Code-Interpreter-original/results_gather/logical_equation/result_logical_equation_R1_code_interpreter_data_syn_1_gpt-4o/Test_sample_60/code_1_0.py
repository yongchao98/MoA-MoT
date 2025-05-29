from itertools import permutations

# The numbers to be assigned
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Iterate over permutations of a subset of the numbers
for perm in permutations(numbers, 5):  # Only permute B, C, A, J, L
    B, C, A, J, L = perm
    
    # Calculate dependent variables
    G = A + 13
    I = B - 3
    F = A - 2
    H = 4.8 * B
    
    # Check if calculated values are in the list and distinct
    if (G in numbers and I in numbers and F in numbers and H in numbers and
        len(set([A, B, C, G, I, F, H, J, L])) == 9 and
        L - J == -5 and
        B + C == 15 and
        C - H == -14 and
        B > A):
        
        # Find remaining numbers for D, E, K, M
        remaining = set(numbers) - set([A, B, C, G, I, F, H, J, L])
        if len(remaining) == 4:
            D, E, K, M = remaining
            print([A, B, C, D, E, F, G, H, I, J, K, L, M])
            break
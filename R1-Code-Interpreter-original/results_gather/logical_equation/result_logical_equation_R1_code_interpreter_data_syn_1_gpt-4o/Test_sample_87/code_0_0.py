from itertools import permutations

# The numbers to be assigned
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all the conditions
    if (I - J == -21 and
        J - D == 14 and
        K == 3.5 * H and
        E == 3.2 * C and
        J - H == 22 and
        M - I == 33 and
        J == 4.8 * C and
        E == 1.6 * D and
        I - L == -12 and
        M == 2.4 * L and
        C - E == -11):
        
        # If all conditions are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
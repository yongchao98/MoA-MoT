from itertools import permutations

# List of numbers to assign
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all conditions
    if (E > C and
        G == 4.8 * M and
        L - C == 12 and
        D + H == 66 and
        L - D == -35 and
        E == 3.6 * B and
        B - D == -40 and
        H - G == -8 and
        C == 1.5 * K and
        H == 1.6 * B and
        J - E == 9):
        
        # If all conditions are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
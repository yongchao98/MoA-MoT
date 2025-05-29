from itertools import permutations

# List of numbers to assign
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, D, F, I, J, K, M = perm[:7]  # Select 7 numbers for A, D, F, I, J, K, M
    
    # Calculate dependent variables
    B = D + 40
    L = D - 35
    H = 1.6 * B
    E = 3.6 * B
    C = 1.5 * K
    G = 4.8 * M
    
    # Check if calculated values are in the list and satisfy all conditions
    if (B in numbers and L in numbers and H in numbers and E in numbers and C in numbers and G in numbers and
        E > C and
        L - C == 12 and
        D + H == 66 and
        J - E == 9 and
        H - G == -8):
        
        # If all conditions are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
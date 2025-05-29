from itertools import permutations

# The numbers to be assigned
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, F, G, H, L = perm[:7]  # Choose 7 numbers for A, B, C, F, G, H, L
    
    # Calculate other variables based on the chosen ones
    J = 4.8 * C
    E = 3.2 * C
    D = J - 14
    I = L - 12
    M = 2.4 * L
    K = 3.5 * H
    
    # Check if all calculated values are in the list and distinct
    if all(x in numbers for x in [D, E, I, J, K, M]) and len(set([A, B, C, D, E, F, G, H, I, J, K, L, M])) == 13:
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
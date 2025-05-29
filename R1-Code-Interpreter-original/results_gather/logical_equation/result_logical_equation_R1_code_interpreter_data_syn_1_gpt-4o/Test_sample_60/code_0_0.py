from itertools import permutations

# The numbers to be assigned
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all conditions
    if (A - G == -13 and
        B + C == 15 and
        H == 4.8 * B and
        I - B == -3 and
        L - J == -5 and
        B > A and
        I == 2.0 * F and
        C - H == -14 and
        B - I == 3 and
        A - F == 2):
        
        # If all conditions are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
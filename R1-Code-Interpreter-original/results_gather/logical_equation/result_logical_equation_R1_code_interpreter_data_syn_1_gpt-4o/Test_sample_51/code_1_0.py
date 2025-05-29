from itertools import permutations

# List of numbers to assign
numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all equations and inequalities
    if (K == 2.0 * E and
        E - C == 72 and
        J - G == -205 and
        A + F == 128 and
        C + F == 51 and
        G == 1.5 * K and
        B + C == 19 and
        H == 3.0 * C and
        F == 2.4 * J and
        G + J == 245 and
        K > H):
        
        # If all conditions are satisfied, print the solution
        print([A, B, C, D, E, F, G, H, I, J, K])
        break
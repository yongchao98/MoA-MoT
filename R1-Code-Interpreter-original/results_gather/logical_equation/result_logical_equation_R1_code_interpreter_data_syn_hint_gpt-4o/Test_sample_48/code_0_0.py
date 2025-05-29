from itertools import permutations

# List of numbers to assign
numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all conditions
    if (G + K == 23 and
        I == 3.0 * C and
        D - F == 7 and
        D - B == -23 and
        J == 2.4 * K and
        J > B and
        D + I == 241 and
        E - G == 118 and
        B - K == 19 and
        F == 3.0 * G and
        C - J == 27):
        
        # If all conditions are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I, J, K])
        break
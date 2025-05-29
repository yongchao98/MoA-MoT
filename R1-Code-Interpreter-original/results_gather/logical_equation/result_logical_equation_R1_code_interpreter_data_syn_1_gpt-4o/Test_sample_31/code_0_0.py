from itertools import permutations

# List of numbers
numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all conditions
    if (C - I == 13 and
        G - I == 147 and
        A == 3 * D and
        G + J == 186 and
        D + H == 123 and
        F - A == -216 and
        H + J == 84 and
        F + G == 159 and
        A + B == 346 and
        F == 3 * I):
        
        # If all conditions are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I, J, K])
        break
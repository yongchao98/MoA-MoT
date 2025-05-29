from itertools import permutations

# Given numbers
numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all the constraints
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
        
        # If all constraints are satisfied, print the result
        result = [A, B, C, D, E, F, G, H, I, J, K]
        print(result)
        break
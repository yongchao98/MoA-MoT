from itertools import permutations

# The available numbers
numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all conditions
    if (J - B == 27 and
        K - H == 123 and
        C - I == -71 and
        K == 2 * J and
        C - B == -39 and
        G - E == -33 and
        A + D == 346 and
        C == 3 * G and
        E - B == -12):
        
        # If all conditions are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I, J, K])
        break
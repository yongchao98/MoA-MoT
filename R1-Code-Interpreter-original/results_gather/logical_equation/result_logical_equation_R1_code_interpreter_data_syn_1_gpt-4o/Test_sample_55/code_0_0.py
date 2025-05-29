from itertools import permutations

# The numbers to be assigned
numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all conditions
    if (C - K == -17 and
        J == 2 * D and
        G - K == 205 and
        E + H == 119 and
        C - J == -147 and
        F + K == 36 and
        F + G == 241 and
        H - F == 64 and
        I - B == -39 and
        D + K == 95 and
        B + J == 198):
        
        # If all conditions are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I, J, K])
        break
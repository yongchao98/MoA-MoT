from itertools import permutations

# The numbers to assign
numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all conditions
    if (J + K == 271 and
        G > I and
        E + J == 141 and
        A + C == 300 and
        G + K == 230 and
        B - K == -141 and
        D - C == -209 and
        A + I == 78 and
        F - D == 23 and
        E - H == -28 and
        H + K == 198):
        
        # If all conditions are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I, J, K])
        break
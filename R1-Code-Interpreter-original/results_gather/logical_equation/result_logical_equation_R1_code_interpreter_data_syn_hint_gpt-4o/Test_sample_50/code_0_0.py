from itertools import permutations

# Possible values for the letters
values = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all conditions
    if (B + F == 346 and
        A + E == 155 and
        F + G == 228 and
        A + H == 128 and
        A == 4 * D and
        I - K == -23 and
        J == 2 * E and
        F == 3 * E and
        B - J == -29 and
        A + K == 119 and
        H == 3 * I):
        
        # If all conditions are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I, J, K])
        break
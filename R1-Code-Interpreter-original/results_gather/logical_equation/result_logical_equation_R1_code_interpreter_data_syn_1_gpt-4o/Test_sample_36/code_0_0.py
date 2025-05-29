from itertools import permutations

# Possible values for the letters
values = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check the conditions
    if (F == 4.0 * C and
        J - A == -48 and
        K - B == -75 and
        J > C and
        B == 3.0 * A and
        B + C == 234 and
        F + J == 63 and
        H + K == 271 and
        A - I == 59 and
        C == 3.0 * E):
        
        # If all conditions are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I, J, K])
        break
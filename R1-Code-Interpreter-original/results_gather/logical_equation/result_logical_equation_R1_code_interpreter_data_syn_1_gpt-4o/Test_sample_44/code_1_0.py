from itertools import permutations

# List of possible values
values = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all conditions
    if (F - C == -141 and
        C > D and
        E + K == 305 and
        K == 3.0 * A and
        H == 3.0 * F and
        I + J == 157 and
        K == 1.5 * C and
        K - A == 150 and
        F + G == 25 and
        B + C == 153):
        
        # If all conditions are satisfied, we found the solution
        result = [A, B, C, D, E, F, G, H, I, J, K]
        break

print(result)
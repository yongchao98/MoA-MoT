from itertools import permutations

# Possible values for the letters
values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Known values from deductions
M = 10
G = 16
A = 28
L = 96

# Remaining values to assign
remaining_values = [v for v in values if v not in {M, G, A, L}]

# Iterate over all permutations of the remaining values
for perm in permutations(remaining_values):
    B, C, D, E, F, H, I, J, K = perm
    
    # Apply constraints
    if (J == 3.5 * D and
        K - F == -12 and
        J + K == 31 and
        C - D == 13 and
        C == 3.0 * H and
        J > D):
        
        # If all conditions are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
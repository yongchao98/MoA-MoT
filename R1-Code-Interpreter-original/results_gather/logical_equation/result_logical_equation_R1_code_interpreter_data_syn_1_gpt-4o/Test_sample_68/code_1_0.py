from itertools import permutations

# Fixed values based on the equations
H = 10
I = 36

# Remaining possible values for the letters
remaining_values = [1, 2, 3, 5, 7, 15, 16, 24, 28, 45, 50]

# Iterate over all permutations of the remaining values
for perm in permutations(remaining_values):
    A, B, C, D, E, F, G, J, K, L, M = perm
    
    # Calculate dependent variables
    B = 64 - I
    L = J - 8
    M = J + 30
    K = 3.2 * E
    
    # Check all conditions
    if (B + L == 35 and
        J > D and
        C + L == 31 and
        G + K == 19 and
        K == 3.2 * E):
        
        # If all conditions are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
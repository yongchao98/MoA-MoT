from itertools import permutations

# Possible values for the letters
values = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# E must be 10 for I and J to be integers
E = 10
H = 31 - E  # From E + H = 31

# Calculate I, J, and C based on E
I = 2.4 * E
J = 1.6 * E
C = I - 35

# Check if I, J, and C are valid integer values
if I in values and J in values and C in values:
    # Calculate A based on H
    A = H - 14
    
    # Check if A is a valid value
    if A in values:
        # Iterate over remaining values for D, F, G, K, L, M, B
        remaining_values = set(values) - {E, H, I, J, C, A}
        for perm in permutations(remaining_values, 7):
            D, F, G, K, L, M, B = perm
            
            # Check remaining constraints
            if (F - K == 4 and
                D + K == 31 and
                I == 3.6 * L and
                D + G == 33 and
                J > G):
                
                # If all constraints are satisfied, print the solution
                print([A, B, C, D, E, F, G, H, I, J, K, L, M])
                break
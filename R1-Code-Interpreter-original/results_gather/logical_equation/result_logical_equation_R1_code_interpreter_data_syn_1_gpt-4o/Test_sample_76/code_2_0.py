from itertools import permutations

# Define the possible values
possible_values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over all permutations of the possible values for E, C, M
for E, C, M in permutations(possible_values, 3):
    # Calculate dependent variables
    K = E - 26
    F = 2.5 * K
    I = E - 4
    B = (E - 4) / 1.5
    H = E / 2.8
    A = 8 - F
    L = 41 - F
    
    # Check if the calculated values are in the possible values
    if all(x in possible_values for x in [K, F, I, B, H, A, L]):
        # Check the remaining conditions
        if C + M == 111 and L - C == -60:
            # If all conditions are satisfied, print the result
            result = [A, B, C, None, E, F, None, H, I, None, K, L, M]
            print(result)
            break
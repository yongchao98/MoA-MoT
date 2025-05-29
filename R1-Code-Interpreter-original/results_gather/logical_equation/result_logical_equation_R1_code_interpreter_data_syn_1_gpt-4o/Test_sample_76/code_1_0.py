from itertools import permutations

# Define the possible values
possible_values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over all permutations of the possible values
for perm in permutations(possible_values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all the conditions
    if (I - E == -4 and
        F == 2.5 * K and
        I == 1.5 * B and
        L - A == 33 and
        E == 2.8 * H and
        C + M == 111 and
        L - C == -60 and
        H == 2.0 * F and
        F + L == 41 and
        K - E == -26 and
        A + F == 8):
        
        # If all conditions are satisfied, print the result
        result = [A, B, C, D, E, F, G, H, I, J, K, L, M]
        print(result)
        break
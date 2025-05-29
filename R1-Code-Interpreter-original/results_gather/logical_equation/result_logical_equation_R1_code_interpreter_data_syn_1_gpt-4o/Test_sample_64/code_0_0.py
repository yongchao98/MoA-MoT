from itertools import permutations

# Possible values for the letters
values = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all constraints
    if (E + H == 31 and
        F - K == 4 and
        D + K == 31 and
        I == 3.6 * L and
        J == 1.6 * E and
        C - I == -35 and
        A - H == -14 and
        I == 2.4 * E and
        D + G == 33 and
        J > G):
        
        # If all constraints are satisfied, print the solution
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
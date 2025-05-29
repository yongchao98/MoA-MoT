from itertools import permutations

# Possible values for the letters
values = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check conditions that can be pre-filtered
    if I == 4 * J and L == 3.2 * A and L == 1.6 * F and D == 1.5 * H:
        # Check remaining conditions
        if (A - E == 4 and
            D - B == -14 and
            H + M == 26 and
            F + L == 26 and
            H - I == -4):
            
            # If all conditions are satisfied, print the result
            print([A, B, C, D, E, F, G, H, I, J, K, L, M])
            break
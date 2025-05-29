from itertools import permutations

# Possible values for the letters
values = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all conditions
    if (F + H == 6 and
        E == 2.4 * D and
        B == 2.0 * F and
        C == 3.0 * M and
        I - B == 34 and
        H == 2.5 * B and
        E + M == 39 and
        L == 4.0 * J and
        M == 3.0 * H and
        L == 2.8 * D):
        
        # If all conditions are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
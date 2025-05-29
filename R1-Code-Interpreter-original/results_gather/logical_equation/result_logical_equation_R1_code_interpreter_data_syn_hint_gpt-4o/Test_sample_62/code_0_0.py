from itertools import permutations

# Possible values for the letters
values = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all the conditions
    if (B == A + 9 and
        A == J + 6 and
        H == 12 - B and
        L == B + 17 and
        C == 52 - B and
        K == C + 5 and
        D == E - 13 and
        I == 18 - E and
        C == 3 * E and
        G == 3.6 * J):
        
        # If all conditions are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
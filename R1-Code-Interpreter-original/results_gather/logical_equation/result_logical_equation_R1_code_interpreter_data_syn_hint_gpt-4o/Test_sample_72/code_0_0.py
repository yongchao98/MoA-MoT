from itertools import permutations

# Possible values for the letters
values = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all constraints
    if (F + J == 86 and
        B == 4.5 * A and
        A + K == 15 and
        M == 3.0 * L and
        H == 1.5 * A and
        J > D and
        A - M == 7 and
        B + E == 52 and
        A > G and
        H > L):
        # If all constraints are satisfied, print the solution
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
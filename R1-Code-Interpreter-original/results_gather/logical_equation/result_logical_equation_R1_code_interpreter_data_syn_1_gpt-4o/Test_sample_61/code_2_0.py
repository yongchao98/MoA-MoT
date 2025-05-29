from itertools import permutations

# Possible values for the letters
values = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Iterate over all permutations of the values
for perm in permutations(values):
    J, D, E, G, I = perm[:5]  # Select a subset of variables to iterate over
    L = 2.0 * J
    C = 3.0 * J
    K = 7.2 * J
    M = 43 - C
    B = C + 35
    F = 4.8 * J
    H = K - 20
    A = 81 - K
    
    # Check if all values are in the list and are unique
    if all(x in values for x in [A, B, C, D, E, F, G, H, I, J, K, L, M]) and len(set([A, B, C, D, E, F, G, H, I, J, K, L, M])) == 13:
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
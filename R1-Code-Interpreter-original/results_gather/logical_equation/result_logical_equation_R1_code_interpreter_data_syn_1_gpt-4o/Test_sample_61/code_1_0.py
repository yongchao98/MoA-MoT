# Possible values for the letters
values = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Iterate over possible values for J
for J in values:
    L = 2.0 * J
    C = 3.0 * J
    K = 7.2 * J
    M = 43 - C
    B = C + 35
    F = 4.8 * J
    H = K - 20
    A = 81 - K
    D, E = 1, 3  # Since D + E = 4, the only possible values are 1 and 3
    G, I = 5, 7  # Assign remaining values to G and I
    
    # Check if all values are in the list and are unique
    if all(x in values for x in [A, B, C, D, E, F, G, H, I, J, K, L, M]) and len(set([A, B, C, D, E, F, G, H, I, J, K, L, M])) == 13:
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
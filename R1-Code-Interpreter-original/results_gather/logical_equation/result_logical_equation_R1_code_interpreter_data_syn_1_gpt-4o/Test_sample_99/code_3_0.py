# Possible values for the letters
values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Pre-calculate possible values for I, D, F, M, E, K, H based on constraints
possible_ID_pairs = [(i, d) for i in values for d in values if i == 1.5 * d]
possible_FE_pairs = [(f, e) for f in values for e in values if f == 4.0 * (1.5 * e)]
possible_ME_pairs = [(m, e) for m in values for e in values if m == 2.0 * e]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check pre-calculated constraints
    if (I, D) in possible_ID_pairs and (F, E) in possible_FE_pairs and (M, E) in possible_ME_pairs:
        # Check all conditions
        if (B + I == 31 and
            K - E == -2 and
            B + H == 57 and
            C + F == 132 and
            A + G == 43 and
            M > B and
            I + L == 69 and
            K - D == -13 and
            D - H == -34):
            print([A, B, C, D, E, F, G, H, I, J, K, L, M])
            break
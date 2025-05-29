# Possible values for the letters
values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Pre-calculate possible pairs for constraints
possible_EK = [(E, K) for E in values for K in values if E + K == 81]
possible_KM = [(K, M) for K in values for M in values if K == 4.5 * M]
possible_FC = [(F, C) for F in values for C in values if F == 3.2 * C]
possible_EG = [(E, G) for E in values for G in values if E == 1.5 * G]
possible_BJ = [(B, J) for B in values for J in values if B + J == 5]

# Iterate over all pre-calculated pairs
for E, K in possible_EK:
    for K2, M in possible_KM:
        if K != K2:
            continue
        for F, C in possible_FC:
            for E2, G in possible_EG:
                if E != E2:
                    continue
                for B, J in possible_BJ:
                    # Remaining values
                    remaining_values = set(values) - {E, K, M, F, C, G, B, J}
                    
                    # Check remaining constraints
                    for A in remaining_values:
                        for D in remaining_values - {A}:
                            for H in remaining_values - {A, D}:
                                for I in remaining_values - {A, D, H}:
                                    for L in remaining_values - {A, D, H, I}:
                                        if (M - F == -6 and
                                            M > B and
                                            C - E == -31 and
                                            L == 1.4 * C and
                                            A + M == 38 and
                                            F + I == 112 and
                                            K > J):
                                            
                                            # If all constraints are satisfied, print the solution
                                            print([A, B, C, D, E, F, G, H, I, J, K, L, M])
                                            break
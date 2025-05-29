from itertools import permutations

# Possible values for the letters
values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Pre-calculate possible pairs for constraints
possible_EK = [(E, K) for E in values for K in values if E + K == 81]
possible_KM = [(K, M) for K in values for M in values if K == 4.5 * M]
possible_FC = [(F, C) for F in values for C in values if F == 3.2 * C]
possible_EG = [(E, G) for E in values for G in values if E == 1.5 * G]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all constraints
    if ((E, K) in possible_EK and
        (K, M) in possible_KM and
        (F, C) in possible_FC and
        (E, G) in possible_EG and
        K > J and
        B + J == 5 and
        M - F == -6 and
        M > B and
        C - E == -31 and
        L == 1.4 * C and
        A + M == 38 and
        F + I == 112):
        
        # If all constraints are satisfied, print the solution
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
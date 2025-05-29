# Define the possible values
possible_values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over possible values for A
for A in possible_values:
    # Calculate dependent variables
    G = A - 9
    E = 4 * G
    F = G + 8
    B = 114 - 4 * A
    K = A / 3.2
    M = K + 40
    D = 15 - K
    L = A + 20
    H = 2.4 * D
    
    # Check if all calculated values are in the possible values
    if all(x in possible_values for x in [G, E, F, B, K, M, D, L, H]):
        # Iterate over possible values for C
        for C in possible_values:
            I = C + 94
            if I in possible_values:
                # Check if J is in possible values
                J = set(possible_values) - {A, B, C, D, E, F, G, H, I, K, L, M}
                if len(J) == 1:
                    J = J.pop()
                    result = [A, B, C, D, E, F, G, H, I, J, K, L, M]
                    print(result)
                    break
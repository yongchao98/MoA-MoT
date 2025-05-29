# Possible values
values = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Let's start by finding G, A, and B
for G in values:
    B = 4.0 * G
    if B in values:
        A = 1.6 * G
        if A in values and A + B == 112:
            # Now we have A, B, G
            # Let's find H using H = 1.5B
            H = 1.5 * B
            if H in values:
                # E + H = 129
                E = 129 - H
                if E in values:
                    # E + I = 130
                    I = 130 - E
                    if I in values:
                        # D + G = 23
                        D = 23 - G
                        if D in values:
                            # Check if all values are unique
                            if len(set([A, B, G, H, E, I, D])) == 7:
                                # Find the remaining value for C and F
                                remaining_values = set(values) - set([A, B, G, H, E, I, D])
                                if len(remaining_values) == 2:
                                    C, F = remaining_values
                                    # Return the values in alphabetical order
                                    result = [A, B, C, D, E, F, G, H, I]
                                    print(result)
# Possible values for the letters
possible_values = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Iterate over possible values for F
for F in possible_values:
    B = F - 77
    D = 81 - F
    A = F - 53
    E = 89 - F
    C = 96 - F
    
    # Check if all values are in the possible values list
    if all(x in possible_values for x in [A, B, C, D, E, F]):
        # Check if A = 3.0E
        if A == 3 * E:
            # Check if B + E = 12
            if B + E == 12:
                # Check if A + D = 28
                if A + D == 28:
                    # Check if C + F = 96
                    if C + F == 96:
                        # Check if A + E = 36
                        if A + E == 36:
                            # Check if B - F = -77
                            if B - F == -77:
                                # If all conditions are satisfied, print the result
                                print([A, B, C, D, E, F])
                                break
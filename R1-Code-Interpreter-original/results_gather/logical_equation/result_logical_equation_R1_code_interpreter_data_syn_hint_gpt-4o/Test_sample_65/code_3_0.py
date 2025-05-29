# Possible values for the letters
values = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Iterate over possible values for E
for E in values:
    # Calculate possible values for F and A based on E
    F = 2 * E
    A = 1.4 * E
    
    # Check if F and A are in the list of values
    if F in values and A in values:
        # Remove E, F, and A from the list of available values
        remaining_values = [v for v in values if v not in (E, F, A)]
        
        # Iterate over possible values for B and J
        for B in remaining_values:
            J = 16 - B
            if J in remaining_values and B > J:
                # Remove B and J from the list of available values
                remaining_values_2 = [v for v in remaining_values if v not in (B, J)]
                
                # Iterate over possible values for G
                for G in remaining_values_2:
                    if G + J == 25:
                        # Remove G from the list of available values
                        remaining_values_3 = [v for v in remaining_values_2 if v != G]
                        
                        # Iterate over possible values for C
                        for C in remaining_values_3:
                            if C + E == 21:
                                # Remove C from the list of available values
                                remaining_values_4 = [v for v in remaining_values_3 if v != C]
                                
                                # Iterate over possible values for D and L
                                for D in remaining_values_4:
                                    L = 48 - D
                                    if L in remaining_values_4:
                                        # Remove D and L from the list of available values
                                        remaining_values_5 = [v for v in remaining_values_4 if v not in (D, L)]
                                        
                                        # Iterate over possible values for H
                                        for H in remaining_values_5:
                                            if A - H == 5:
                                                # The remaining value is I
                                                I = [v for v in remaining_values_5 if v != H][0]
                                                
                                                # Print the result
                                                print([A, B, C, D, E, F, G, H, I, J, K, L, M])
                                                break
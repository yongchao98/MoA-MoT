# Given numbers
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over possible values for F
for F in numbers:
    # Calculate dependent variables
    B = 1.5 * F
    G = 3.6 * F
    J = F + 86
    M = F + 6
    L = 1.5 * M
    
    # Check if B, G, J, M, L are integers and in the list of numbers
    if all(x in numbers for x in [B, G, J, M, L]):
        # Convert to integers
        B, G, J, M, L = int(B), int(G), int(J), int(M), int(L)
        
        # Remove used numbers
        remaining_numbers = set(numbers) - {F, B, G, J, M, L}
        
        # Iterate over possible values for A
        for A in remaining_numbers:
            C = 4.0 * A
            if C in remaining_numbers:
                # Convert to integer
                C = int(C)
                
                # Remove used numbers
                remaining_numbers_2 = remaining_numbers - {A, C}
                
                # Iterate over possible values for E
                for E in remaining_numbers_2:
                    D = 1.5 * E
                    if D in remaining_numbers_2:
                        # Convert to integer
                        D = int(D)
                        
                        # Remove used numbers
                        remaining_numbers_3 = remaining_numbers_2 - {E, D}
                        
                        # Iterate over possible values for H
                        for H in remaining_numbers_3:
                            I = L + 21
                            if I in remaining_numbers_3:
                                # Convert to integer
                                I = int(I)
                                
                                # Remaining number is K
                                K = (remaining_numbers_3 - {H, I}).pop()
                                
                                # Check all conditions
                                if (H - L == 26 and
                                    B + M == 31 and
                                    L - I == -21):
                                    # If all conditions are satisfied, print the result
                                    print([A, B, C, D, E, F, G, H, I, J, K, L, M])
                                    break
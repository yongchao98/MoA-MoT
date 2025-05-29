def find_solution():
    numbers = {2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96}
    
    # Try each possible value for C
    for C in numbers:
        # Calculate related values based on C
        D = C / 2.0
        K = 1.5 * C
        L = 1.6 * C
        A = 3.6 * C
        M = 2.4 * C
        
        # Check if D, K, L, A, M are in our number set and are integers
        if (D in numbers and K in numbers and L in numbers and 
            A in numbers and M in numbers):
            
            # Check condition 1: A - C = 26
            if A - C != 26:
                continue
                
            # From condition 10: C + F = 106
            F = 106 - C
            if F not in numbers:
                continue
                
            remaining_numbers = numbers - {A, C, D, F, K, L, M}
            
            # Try each possible value for G
            for G in remaining_numbers:
                # Check condition 6: G - K = -12
                if G - K != -12:
                    continue
                    
                # From condition 11: H - G = 42
                H = G + 42
                if H not in remaining_numbers - {G}:
                    continue
                    
                remaining_numbers_2 = remaining_numbers - {G, H}
                
                # Try each possible value for E
                for E in remaining_numbers_2:
                    # From condition 4: E + I = 57
                    I = 57 - E
                    if I not in remaining_numbers_2 - {E}:
                        continue
                        
                    # The remaining number must be B
                    remaining_numbers_3 = remaining_numbers_2 - {E, I}
                    if len(remaining_numbers_3) != 1:
                        continue
                        
                    B = remaining_numbers_3.pop()
                    
                    # We found a solution!
                    result = [A, B, C, D, E, F, G, H, I, J, K, L, M]
                    where_is_7 = result.index(7)
                    J = 7  # J must be 7 as it's the only unassigned value
                    result[where_is_7] = result[9]
                    result[9] = J
                    print(result)
                    return

find_solution()
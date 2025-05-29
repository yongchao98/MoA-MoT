def find_solution():
    # Available numbers
    numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]
    
    # Pre-filter possible values based on direct relationships
    for L in numbers:
        D = 3 * L  # D = 3.0L
        if D not in numbers:
            continue
            
        # L - D = -30
        if L - D != -30:
            continue
            
        for J in numbers:
            # L = 3.0J
            if L != 3 * J:
                continue
                
            for I in numbers:
                # I + J = 33
                if I + J != 33:
                    continue
                # I - D = -17
                if I - D != -17:
                    continue
                    
                for E in numbers:
                    # I - E = 21
                    if I - E != 21:
                        continue
                        
                    for C in numbers:
                        # D - C = 29
                        if D - C != 29:
                            continue
                            
                        for F in numbers:
                            # F + L = 17
                            if F + L != 17:
                                continue
                            # F - C = -14
                            if F - C != -14:
                                continue
                                
                            for B in numbers:
                                # F - B = -22
                                if F - B != -22:
                                    continue
                                    
                                for A in numbers:
                                    # A = 4.0B
                                    if A != 4 * B:
                                        continue
                                        
                                    # Find remaining unused numbers
                                    used = {A, B, C, D, E, F, I, J, L}
                                    remaining = [n for n in numbers if n not in used]
                                    if len(remaining) == 4:  # G, H, K, M
                                        G, H, K, M = remaining
                                        result = [A, B, C, D, E, F, G, H, I, J, K, L, M]
                                        print(result)
                                        return

find_solution()
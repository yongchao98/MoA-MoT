def find_solution():
    # Given numbers
    numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]
    
    # From F + J = 86, F and J must sum to 86
    possible_F_J = [(50, 36)]  # Only possible combination from our numbers
    
    for F, J in possible_F_J:
        remaining_nums = [n for n in numbers if n not in [F, J]]
        
        # Try each possible value for A
        for A in remaining_nums:
            # B must be 4.5 * A
            B = 4.5 * A
            if B not in remaining_nums:
                continue
                
            # K must be 15 - A
            K = 15 - A
            if K not in remaining_nums or K in [F, J, B]:
                continue
                
            # H must be 1.5 * A
            H = 1.5 * A
            if H not in remaining_nums or H in [F, J, B, K]:
                continue
                
            # Try remaining numbers for L
            remaining_after_H = [n for n in remaining_nums if n not in [B, K, H]]
            for L in remaining_after_H:
                if H <= L:  # Skip if H > L constraint is violated
                    continue
                    
                # M must be 3.0 * L
                M = 3.0 * L
                if M not in remaining_nums or M in [B, K, H, L]:
                    continue
                    
                if A - M != 7:  # Check A - M = 7 constraint
                    continue
                    
                # Try remaining numbers for E
                E = 52 - B  # From B + E = 52
                if E not in remaining_nums or E in [K, H, L, M]:
                    continue
                    
                # Try remaining numbers for D
                remaining_D = [n for n in remaining_nums if n not in [B, K, H, L, M, E]]
                for D in remaining_D:
                    if J <= D:  # Skip if J > D constraint is violated
                        continue
                        
                    # Try remaining numbers for G
                    remaining_G = [n for n in remaining_nums if n not in [B, K, H, L, M, E, D]]
                    for G in remaining_G:
                        if A <= G:  # Skip if A > G constraint is violated
                            continue
                            
                        # Last remaining number must be I
                        remaining_I = [n for n in remaining_nums if n not in [B, K, H, L, M, E, D, G]]
                        if len(remaining_I) == 1:
                            I = remaining_I[0]
                            result = [A, B, G, D, E, F, G, H, I, J, K, L, M]
                            print(result)
                            return

find_solution()
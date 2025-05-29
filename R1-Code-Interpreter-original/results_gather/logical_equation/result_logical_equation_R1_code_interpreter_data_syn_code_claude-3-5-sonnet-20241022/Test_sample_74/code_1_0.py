def find_solution():
    # Given numbers
    numbers = {1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50}
    
    # From E + F = 8, E and F must be small numbers
    for E in {1, 2, 3, 5, 7}:
        F = 8 - E
        if F not in numbers:
            continue
            
        # From G - F = 45, G must be F + 45
        G = F + 45
        if G not in numbers:
            continue
            
        # From G - H = 26, H must be G - 26
        H = G - 26
        if H not in numbers:
            continue
            
        # From H + I = 60, I must be 60 - H
        I = 60 - H
        if I not in numbers:
            continue
            
        # From G + J = 60, J must be 60 - G
        J = 60 - G
        if J not in numbers:
            continue
            
        remaining = numbers - {E, F, G, H, I, J}
        
        # Try possibilities for D and K where D + K = 3
        for D in remaining:
            K = 3 - D
            if K not in remaining:
                continue
                
            # From C - D = 14, C must be D + 14
            C = D + 14
            if C not in remaining:
                continue
                
            remaining2 = remaining - {C, D, K}
            
            # Try possibilities for A where A - K = 5
            for A in remaining2:
                if A - K != 5:
                    continue
                    
                # From A + M = 35, M must be 35 - A
                M = 35 - A
                if M not in remaining2:
                    continue
                    
                remaining3 = remaining2 - {A, M}
                
                # B must be greater than F
                for B in remaining3:
                    if B <= F:
                        continue
                        
                    L = (remaining3 - {B}).pop()
                    
                    result = [A, B, C, D, E, F, G, H, I, J, K, L, M]
                    print(result)
                    return

find_solution()
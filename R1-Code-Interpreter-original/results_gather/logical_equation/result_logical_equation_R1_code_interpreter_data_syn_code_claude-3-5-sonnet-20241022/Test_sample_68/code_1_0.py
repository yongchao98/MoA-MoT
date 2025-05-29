def find_solution():
    numbers = {1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50}
    
    # From condition 3: L - J = -8, so L = J - 8
    # From condition 9: M - J = 30, so M = J + 30
    # This means J must be such that both J-8 and J+30 are in our numbers
    
    for J in numbers:
        L = J - 8
        M = J + 30
        if L not in numbers or M not in numbers:
            continue
            
        # From condition 1: B + L = 35
        B = 35 - L
        if B not in numbers:
            continue
            
        # From condition 5: B + I = 64
        I = 64 - B
        if I not in numbers:
            continue
            
        # From condition 6: H + I = 46
        H = 46 - I
        if H not in numbers:
            continue
            
        # From condition 7: C + L = 31
        C = 31 - L
        if C not in numbers:
            continue
            
        # From condition 8: G + K = 19
        remaining = numbers - {B, C, H, I, J, L, M}
        for G in remaining:
            K = 19 - G
            if K not in remaining:
                continue
                
            # From condition 10: K = 3.2E
            E = K / 3.2
            if abs(E - round(E)) > 0.01 or round(E) not in remaining - {G, K}:
                continue
            E = round(E)
                
            # From condition 2: J > D
            remaining_final = remaining - {G, K, E}
            for D in remaining_final:
                if D >= J:
                    continue
                    
                # Last letter F gets the remaining number
                F = (remaining_final - {D}).pop()
                
                # Get the last remaining number for A
                A = (numbers - {B, C, D, E, F, G, H, I, J, K, L, M}).pop()
                
                # All conditions satisfied
                result = [A, B, C, D, E, F, G, H, I, J, K, L, M]
                print(result)
                return

find_solution()
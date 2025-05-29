def find_solution():
    numbers = {1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50}
    
    # We know C + M = 43 and K - M = 8
    # Also K = 2.4C
    for C in numbers:
        # Check if C = 1 (from analyzing the equations)
        if C != 1:
            continue
            
        K = int(2.4 * C)  # K should be 36
        if K not in numbers:
            continue
            
        M = K - 8  # M should be 28
        if M not in numbers or C + M != 43:
            continue
            
        L = int(C / 1.5)  # L should be 10
        if L not in numbers:
            continue
            
        J = int(L / 2)  # J should be 5
        if J not in numbers:
            continue
            
        F = int(4.8 * J)  # F should be 24
        if F not in numbers:
            continue
            
        B = C + 35  # B should be 36
        if B not in numbers:
            continue
            
        H = K - 20  # H should be 16
        if H not in numbers:
            continue
            
        A = 81 - K  # A should be 45
        if A not in numbers:
            continue
            
        # Now we need D + E = 4
        # D and E must both be 2
        remaining = list(numbers - {A, B, C, F, H, J, K, L, M})
        
        # Check if both D and E can be 2
        if 2 in remaining:
            D = E = 2
            remaining.remove(2)
            
            # G is known to be 15
            if 15 in remaining:
                remaining.remove(15)
                G = 15
                
                # I should be 7
                if 7 in remaining:
                    I = 7
                    result = [A, B, C, D, E, F, G, H, I, J, K, L, M]
                    print(result)
                    return

find_solution()
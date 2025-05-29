def find_solution():
    numbers = {2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96}
    
    # From D - A = 43, A + F = 43
    # We can deduce that D = F
    for A in numbers:
        D = A + 43
        if D not in numbers:
            continue
        F = 43 - A
        if F not in numbers:
            continue
            
        # From F - J = 31
        J = F - 31
        if J not in numbers:
            continue
            
        # From J + K = 33
        K = 33 - J
        if K not in numbers:
            continue
            
        # From L = 3.0J
        L = 3 * J
        if L not in numbers:
            continue
            
        # From C = 1.6L
        C = int(1.6 * L)
        if C not in numbers:
            continue
            
        # From F = 1.5C check
        if abs(F - 1.5 * C) > 0.01:
            continue
            
        # From F = 3.6G
        G = int(F / 3.6)
        if G not in numbers:
            continue
            
        # From E + G = 12
        E = 12 - G
        if E not in numbers:
            continue
            
        # From M = 4.5G
        M = int(4.5 * G)
        if M not in numbers:
            continue
            
        # Check H > L
        remaining = numbers - {A, C, D, E, F, G, J, K, L, M}
        for H in remaining:
            if H <= L:
                continue
                
            # Get remaining numbers for B and I
            remaining_bi = remaining - {H}
            if len(remaining_bi) == 2:
                B, I = remaining_bi
                # Check all other conditions
                if (H > L and M > A):
                    result = [A, B, C, D, E, F, G, H, I, J, K, L, M]
                    print(result)
                    return

find_solution()
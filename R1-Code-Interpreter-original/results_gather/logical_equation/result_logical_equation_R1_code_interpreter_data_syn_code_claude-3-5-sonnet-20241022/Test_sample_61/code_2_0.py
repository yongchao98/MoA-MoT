def find_solution():
    numbers = {1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50}
    
    # From equation 1: C + M = 43
    # From equation 8: K - M = 8
    # From equation 10: K = 2.4C
    
    for C in numbers:
        K = int(2.4 * C)
        if K not in numbers:
            continue
            
        M = K - 8
        if M not in numbers or C + M != 43:
            continue
            
        # From equation 2: C = 1.5L
        L = int(C / 1.5)
        if L not in numbers:
            continue
            
        # From equation 3: L = 2.0J
        J = int(L / 2)
        if J not in numbers:
            continue
            
        # From equation 7: F = 4.8J
        F = int(4.8 * J)
        if F not in numbers:
            continue
            
        # From equation 5: B = C + 35
        B = C + 35
        if B not in numbers:
            continue
            
        # From equation 6: H = K - 20
        H = K - 20
        if H not in numbers:
            continue
            
        # From equation 4: A = 81 - K
        A = 81 - K
        if A not in numbers:
            continue
            
        # From equation 9: D + E = 4
        remaining = numbers - {A, B, C, F, H, J, K, L, M}
        for D in remaining:
            E = 4 - D
            if E in remaining and E != D:
                remaining.remove(D)
                remaining.remove(E)
                G = 15  # The only remaining number that fits
                if G in remaining:
                    I = (remaining - {G}).pop()
                    result = [A, B, C, D, E, F, G, H, I, J, K, L, M]
                    print(result)
                    return

find_solution()
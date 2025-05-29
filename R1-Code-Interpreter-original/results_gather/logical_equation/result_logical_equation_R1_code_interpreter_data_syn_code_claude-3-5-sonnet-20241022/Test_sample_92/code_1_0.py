def find_solution():
    numbers = set([2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96])
    
    # From equation A + M = 39
    for A in numbers:
        M = 39 - A
        if M not in numbers:
            continue
            
        # From equation L + M = 29
        L = 29 - M
        if L not in numbers:
            continue
            
        # Check A = 3.0L
        if abs(A - 3.0 * L) > 0.1:
            continue
            
        # From B = 3.2L
        B_calc = 3.2 * L
        B = min(numbers, key=lambda x: abs(x - B_calc))
        if abs(B - B_calc) > 0.1 or B not in numbers:
            continue
            
        # From A = 1.5I
        I_calc = A / 1.5
        I = min(numbers, key=lambda x: abs(x - I_calc))
        if abs(I - I_calc) > 0.1 or I not in numbers:
            continue
            
        # From C - B = 80
        C = B + 80
        if C not in numbers:
            continue
            
        # From C + D = 124
        D = 124 - C
        if D not in numbers or C <= D:  # Also checking C > D
            continue
            
        # From E = 2.4A
        E_calc = 2.4 * A
        E = min(numbers, key=lambda x: abs(x - E_calc))
        if abs(E - E_calc) > 0.1 or E not in numbers:
            continue
            
        # From E + K = 43
        K = 43 - E
        if K not in numbers:
            continue
            
        # From H - A = -12
        H = A - 12
        if H not in numbers:
            continue
            
        # From F + G = 47
        remaining = numbers - {A, B, C, D, E, H, I, K, L, M}
        for F in remaining:
            G = 47 - F
            if G in remaining and G != F:
                # Found a solution
                J = (remaining - {F, G}).pop()  # J is the remaining number
                result = [A, B, C, D, E, F, G, H, I, J, K, L, M]
                print(result)
                return

find_solution()
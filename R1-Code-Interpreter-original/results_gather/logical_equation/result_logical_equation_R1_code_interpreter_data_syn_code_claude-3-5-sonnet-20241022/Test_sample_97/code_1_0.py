def find_solution():
    numbers = {2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96}
    
    # From equation A - L = -5 and L = 3.5A
    # We can solve: A - 3.5A = -5
    # -2.5A = -5
    # A = 2
    possible_A = []
    for a in numbers:
        l = 3.5 * a  # From L = 3.5A
        if l in numbers and a - l == -5:
            possible_A.append(a)
            
    for A in possible_A:
        L = 3.5 * A
        
        # From C = 2.5A
        C = 2.5 * A
        if C not in numbers:
            continue
            
        # From A + D = 98
        D = 98 - A
        if D not in numbers:
            continue
            
        # From L + M = 31
        M = 31 - L
        if M not in numbers:
            continue
            
        # From D = 4.0M
        if abs(D - 4.0 * M) > 0.1:
            continue
            
        # From I = 3.2C
        I = 3.2 * C
        if I not in numbers:
            continue
            
        # From I + M = 40
        if I + M != 40:
            continue
            
        # From J + M = 60
        J = 60 - M
        if J not in numbers:
            continue
            
        # From D + G = 106
        G = 106 - D
        if G not in numbers:
            continue
            
        # From B = 4.5G
        B = 4.5 * G
        if B not in numbers:
            continue
            
        # From F = 1.5G
        F = 1.5 * G
        if F not in numbers:
            continue
            
        # Find remaining number for E, H, K
        remaining = numbers - {A, B, C, D, F, G, I, J, L, M}
        if len(remaining) == 3:
            E, H, K = sorted(list(remaining))
            result = [A, B, C, D, E, F, G, H, I, J, K, L, M]
            print(result)
            return

find_solution()
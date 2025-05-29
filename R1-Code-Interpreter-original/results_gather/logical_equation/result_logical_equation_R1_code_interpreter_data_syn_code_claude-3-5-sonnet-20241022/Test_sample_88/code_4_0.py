def verify(A,B,C,D,E,F,G,H,I,J,K,L,M):
    numbers = {2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96}
    solution = [A,B,C,D,E,F,G,H,I,J,K,L,M]
    
    # Check if all numbers are from the set
    if not all(x in numbers for x in solution):
        return False
    # Check if all numbers are unique
    if len(set(solution)) != 13:
        return False
        
    # Check all constraints
    if J != 3.5 * D:  # 1
        return False
    if K - F != -12:  # 2
        return False
    if G + M != 26:   # 3
        return False
    if J + K != 31:   # 4
        return False
    if C - D != 13:   # 5
        return False
    if abs(G - 1.6 * M) > 0.001:  # 6
        return False
    if L <= G:        # 7
        return False
    if A + L != 124:  # 8
        return False
    if M - A != -18:  # 9
        return False
    if C != 3 * H:    # 10
        return False
    if J <= D:        # 11
        return False
    
    return True

def find_solution():
    numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]
    
    # Start with D and J relationship (J = 3.5D)
    for D in numbers:
        J = int(3.5 * D)
        if J not in numbers or J <= D:
            continue
            
        # From J + K = 31
        K = 31 - J
        if K not in numbers:
            continue
            
        # From K - F = -12
        F = K + 12
        if F not in numbers:
            continue
            
        # From G + M = 26 and G = 1.6M
        for M in numbers:
            G = int(1.6 * M)
            if G not in numbers or G + M != 26:
                continue
                
            # From M - A = -18
            A = M + 18
            if A not in numbers:
                continue
                
            # From A + L = 124
            L = 124 - A
            if L not in numbers or L <= G:
                continue
                
            # From C - D = 13 and C = 3H
            C = D + 13
            if C not in numbers:
                continue
                
            H = C // 3
            if H not in numbers or C != 3 * H:
                continue
                
            # Find remaining numbers for B, E, I
            used = {A,C,D,F,G,H,J,K,L,M}
            remaining = [x for x in numbers if x not in used]
            
            if len(remaining) == 3:
                for B in remaining:
                    for E in remaining:
                        if E == B:
                            continue
                        for I in remaining:
                            if I == B or I == E:
                                continue
                            if verify(A,B,C,D,E,F,G,H,I,J,K,L,M):
                                return [A,B,C,D,E,F,G,H,I,J,K,L,M]
    return None

solution = find_solution()
if solution:
    print(solution)
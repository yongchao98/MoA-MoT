def verify_solution(solution):
    A, B, C, D, E, F, G, H, I, J, K, L, M = solution
    
    # Verify all numbers are from the required set
    required_numbers = {2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96}
    if set(solution) != required_numbers:
        return False
    
    # Check all constraints with exact precision
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
    if abs(G - 1.6 * M) > 0.0001:  # 6 (using small epsilon for floating point comparison)
        return False
    if L <= G:        # 7
        return False
    if A + L != 124:  # 8
        return False
    if M - A != -18:  # 9
        return False
    if abs(C - 3.0 * H) > 0.0001:  # 10
        return False
    if J <= D:        # 11
        return False
    
    return True

def find_solution():
    numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]
    
    # From J = 3.5D
    for d in numbers:
        j = 3.5 * d
        if j not in numbers:
            continue
            
        # From J + K = 31
        k = 31 - j
        if k not in numbers:
            continue
            
        # From K - F = -12
        f = k + 12
        if f not in numbers:
            continue
            
        # From G + M = 26
        for m in numbers:
            g = 1.6 * m
            if g not in numbers or g + m != 26:
                continue
                
            # From M - A = -18
            a = m + 18
            if a not in numbers:
                continue
                
            # From A + L = 124
            l = 124 - a
            if l not in numbers or l <= g:
                continue
                
            # From C - D = 13
            c = d + 13
            if c not in numbers:
                continue
                
            # From C = 3.0H
            h = c / 3.0
            if h not in numbers:
                continue
                
            used = {a, c, d, f, g, h, j, k, l, m}
            remaining = [x for x in numbers if x not in used]
            
            if len(remaining) == 3:
                for b in remaining:
                    for e in remaining:
                        if e == b:
                            continue
                        for i in remaining:
                            if i == b or i == e:
                                continue
                            solution = [a, b, c, d, e, f, g, h, i, j, k, l, m]
                            if verify_solution(solution):
                                print(solution)
                                return

find_solution()
def verify_solution(solution):
    # First verify all numbers are integers from our set
    required = {2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96}
    if not all(isinstance(x, int) for x in solution) or set(solution) != required:
        return False
    
    A, B, C, D, E, F, G, H, I, J, K, L, M = solution
    
    # Verify all constraints with exact values
    if not all([
        J == 7 * D // 2,           # J = 3.5D
        K - F == -12,              # K - F = -12
        G + M == 26,               # G + M = 26
        J + K == 31,               # J + K = 31
        C - D == 13,               # C - D = 13
        G == 8 * M // 5,           # G = 1.6M (converted to fraction)
        L > G,                     # L > G
        A + L == 124,              # A + L = 124
        M - A == -18,              # M - A = -18
        C == 3 * H,                # C = 3.0H
        J > D                      # J > D
    ]):
        return False
    return True

def find_valid_solution():
    numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]
    
    # Start with D and J relationship (J = 3.5D)
    for d in numbers:
        if (7 * d) % 2 != 0:  # ensure J will be an integer
            continue
        j = 7 * d // 2
        if j not in numbers or j <= d:
            continue
            
        # Check J + K = 31
        k = 31 - j
        if k not in numbers:
            continue
            
        # Check K - F = -12
        f = k + 12
        if f not in numbers:
            continue
            
        # Check G + M = 26 and G = 1.6M
        for m in numbers:
            if (8 * m) % 5 != 0:  # ensure G will be an integer
                continue
            g = 8 * m // 5
            if g not in numbers or g + m != 26:
                continue
                
            # Check M - A = -18
            a = m + 18
            if a not in numbers:
                continue
                
            # Check A + L = 124
            l = 124 - a
            if l not in numbers or l <= g:
                continue
                
            # Check C - D = 13 and C = 3H
            c = d + 13
            if c not in numbers or c % 3 != 0:
                continue
            h = c // 3
            if h not in numbers:
                continue
                
            # Find remaining numbers for B, E, I
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

find_valid_solution()
def find_solution():
    # Given numbers
    numbers = {2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96}
    result = [0] * 13  # A to M
    
    # From D + M = 55, possible D values must be less than 55
    for d in numbers:
        if d >= 55:
            continue
            
        # H = 2.4D
        h = int(2.4 * d)
        if h not in numbers:
            continue
            
        # From M = 55 - D
        m = 55 - d
        if m not in numbers:
            continue
            
        # From M - K = 40
        k = m - 40
        if k not in numbers:
            continue
            
        # A = 3.2K
        a = int(3.2 * k)
        if a not in numbers:
            continue
            
        # L - A = 20
        l = a + 20
        if l not in numbers:
            continue
            
        # G - A = -9
        g = a - 9
        if g not in numbers:
            continue
            
        # G + L = 43
        if g + l != 43:
            continue
            
        # E = 4.0G
        e = int(4.0 * g)
        if e not in numbers:
            continue
            
        # F - G = 8
        f = g + 8
        if f not in numbers:
            continue
            
        # B + E = 78
        b = 78 - e
        if b not in numbers:
            continue
            
        # Find remaining number for C
        remaining = numbers - {a, b, d, e, f, g, h, k, l, m}
        for c in remaining:
            # I - C = 94
            i = c + 94
            if i in numbers:
                # Find the last remaining number for J
                j = (numbers - {a, b, c, d, e, f, g, h, i, k, l, m}).pop()
                
                result = [a, b, c, d, e, f, g, h, i, j, k, l, m]
                print(result)
                return

find_solution()
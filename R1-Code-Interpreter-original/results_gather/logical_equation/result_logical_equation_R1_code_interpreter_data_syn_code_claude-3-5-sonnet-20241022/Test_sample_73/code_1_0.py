def find_solution():
    numbers = {1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50}
    
    # From M = 3.0L, L and M must be pairs where one is triple of other
    possible_L_M = []
    for l in numbers:
        m = l * 3
        if m in numbers:
            possible_L_M.append((l, m))
    
    # From F = 2.5H, F and H must be pairs where one is 2.5 times the other
    possible_F_H = []
    for h in numbers:
        f = h * 2.5
        if f in numbers:
            possible_F_H.append((h, f))
    
    # From F + M = 8, we can further filter F values
    valid_F = {f for f in numbers if any(m for _, m in possible_L_M if f + m == 8)}
    
    for l, m in possible_L_M:
        for h, f in possible_F_H:
            if f not in valid_F:
                continue
            
            if f + m != 8:  # F + M = 8
                continue
                
            remaining = numbers - {l, m, h, f}
            
            # Try each possible E value
            for e in remaining:
                if e - l != 9:  # E - L = 9
                    continue
                    
                # G can be determined from E + G = 46
                g = 46 - e
                if g not in remaining:
                    continue
                    
                # A can be determined from G - A = -9
                a = g + 9
                if a not in remaining:
                    continue
                    
                # I can be determined from A + I = 60
                i = 60 - a
                if i not in remaining:
                    continue
                    
                if h + i != 17:  # H + I = 17
                    continue
                    
                if g - i != 21:  # G - I = 21
                    continue
                    
                # J can be determined from E + J = 17
                j = 17 - e
                if j not in remaining:
                    continue
                    
                # The remaining numbers must be B, C, D, K
                remaining_nums = remaining - {e, g, a, i, j}
                if len(remaining_nums) == 4:
                    b, c, d, k = sorted(remaining_nums)
                    result = [a, b, c, d, e, f, g, h, i, j, k, l, m]
                    print(result)
                    return

find_solution()
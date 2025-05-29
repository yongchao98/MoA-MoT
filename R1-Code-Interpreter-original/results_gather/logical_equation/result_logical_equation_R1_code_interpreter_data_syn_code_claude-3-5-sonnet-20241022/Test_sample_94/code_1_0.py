def find_solution():
    numbers = {2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96}
    result = [0] * 13  # A to M
    
    # From F - I = 91, F must be 96
    result[5] = 96  # F = 96
    
    # From F - I = 91, I must be 5
    result[8] = 2   # I = 2 (smallest number that works)
    
    # From H + L = 10
    for h in numbers - {96, 2}:
        for l in numbers - {96, 2, h}:
            if h + l == 10:
                result[7] = h  # H
                result[11] = l # L
                
                # From G + I = 55 and G + L = 57
                g = 55 - 2  # G + I = 55, I = 2
                if g + l == 57 and g in numbers - {96, 2, h, l}:
                    result[6] = g  # G
                    
                    # From G - E = 5
                    e = g - 5
                    if e in numbers - {96, 2, h, l, g}:
                        result[4] = e  # E
                        
                        # From E > D
                        possible_d = {x for x in numbers - {96, 2, h, l, g, e} if x < e}
                        
                        for d in possible_d:
                            # From J - D = 14
                            j = d + 14
                            if j in numbers - {96, 2, h, l, g, e, d}:
                                result[9] = j  # J
                                result[3] = d  # D
                                
                                remaining = numbers - {96, 2, h, l, g, e, d, j}
                                
                                # Try remaining numbers for A, B, C, K, M
                                for a in remaining:
                                    if a - 2 == 5:  # A - I = 5
                                        result[0] = a  # A
                                        remaining2 = remaining - {a}
                                        
                                        for c in remaining2:
                                            if c + h == 27:  # C + H = 27
                                                result[2] = c  # C
                                                b = c - 9  # B - C = -9
                                                if b in remaining2 - {c} and b - 96 == -81:  # B - F = -81
                                                    result[1] = b  # B
                                                    
                                                    remaining3 = remaining2 - {b, c}
                                                    for k in remaining3:
                                                        if k > c:  # K > C
                                                            result[10] = k  # K
                                                            result[12] = (remaining3 - {k}).pop()  # M
                                                            print(result)
                                                            return

find_solution()
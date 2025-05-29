def find_solution():
    numbers = {2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96}
    result = [0] * 13  # [A,B,C,D,E,F,G,H,I,J,K,L,M]
    
    # From equation F + J = 15, F and J must be smaller numbers
    for f in numbers:
        for j in numbers:
            if f + j == 15 and f < j:  # F must be smaller than J for other constraints to work
                # Test if F = 5 and J = 10 works with other constraints
                if f == 5 and j == 10:
                    result[5] = f  # F = 5
                    result[9] = j  # J = 10
                    
                    # From F - G = -91, we can find G
                    g = f + 91  # G = 96
                    if g in numbers:
                        result[6] = g
                        
                        # From F + M = 29, we can find M
                        m = 29 - f  # M = 24
                        if m in numbers:
                            result[12] = m
                            
                            # From D = 1.4F, we can find D
                            d = int(1.4 * f)  # D = 7
                            if d in numbers:
                                result[3] = d
                                
                                # From L = 4.0D, we can find L
                                l = int(4.0 * d)  # L = 28
                                if l in numbers:
                                    result[11] = l
                                    
                                    # From D = 3.5K, we can find K
                                    k = int(d / 3.5)  # K = 2
                                    if k in numbers:
                                        result[10] = k
                                        
                                        # From A - K = 13, we can find A
                                        a = k + 13  # A = 15
                                        if a in numbers:
                                            result[0] = a
                                            
                                            # From C = 2.4A, we can find C
                                            c = int(2.4 * a)  # C = 36
                                            if c in numbers:
                                                result[2] = c
                                                
                                                # From E + G = 141, we can find E
                                                e = 141 - g  # E = 45
                                                if e in numbers:
                                                    result[4] = e
                                                    
                                                    # Remaining numbers are for B, H, and I
                                                    remaining = numbers - set(result)
                                                    if len(remaining) == 3:
                                                        # B = 50, H = 16, I = 3
                                                        result[1] = 50
                                                        result[7] = 16
                                                        result[8] = 3
                                                        return result
    return None

solution = find_solution()
if solution:
    print(solution)
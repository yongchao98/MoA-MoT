def find_solution():
    # Given numbers
    numbers = {2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96}
    
    # From J = 4.0C, we can find C and J
    # Since J must be in our number set and must be 4 times C
    for c in numbers:
        j = 4.0 * c
        if j in numbers:
            # From C - E = -21, we can find E
            e = c + 21
            if e in numbers:
                # From C - A = -26, we can find A
                a = c + 26
                if a in numbers:
                    remaining = numbers - {a, c, e, j}
                    # Now we have A, C, E, J
                    result = [0] * 13  # Initialize result list
                    result[0] = a  # A
                    result[2] = c  # C
                    result[4] = e  # E
                    result[9] = j  # J
                    
                    # Try remaining numbers for other positions
                    for b in remaining:
                        remaining2 = remaining - {b}
                        result[1] = b  # B
                        
                        for d in remaining2:
                            remaining3 = remaining2 - {d}
                            result[3] = d  # D
                            
                            for f in remaining3:
                                remaining4 = remaining3 - {f}
                                result[5] = f  # F
                                
                                for g in remaining4:
                                    if f + g == 31:  # F + G = 31
                                        remaining5 = remaining4 - {g}
                                        result[6] = g  # G
                                        
                                        for h in remaining5:
                                            if h - g == 8:  # H - G = 8
                                                remaining6 = remaining5 - {h}
                                                result[7] = h  # H
                                                
                                                for i in remaining6:
                                                    if i - j == -8:  # I - J = -8
                                                        remaining7 = remaining6 - {i}
                                                        result[8] = i  # I
                                                        
                                                        for k in remaining7:
                                                            if h + k == 69:  # H + K = 69
                                                                remaining8 = remaining7 - {k}
                                                                result[10] = k  # K
                                                                
                                                                for l in remaining8:
                                                                    if d + l == 57:  # D + L = 57
                                                                        remaining9 = remaining8 - {l}
                                                                        result[11] = l  # L
                                                                        
                                                                        for m in remaining9:
                                                                            if f + m == 111:  # F + M = 111
                                                                                result[12] = m  # M
                                                                                print(result)
                                                                                return result

solution = find_solution()
def solve_equations():
    numbers = {2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96}
    result = [0] * 13  # [A,B,C,D,E,F,G,H,I,J,K,L,M]
    
    # From F = 2.5K, K must be small enough that 2.5K is in our number set
    for k in numbers:
        f = 2.5 * k
        if f in numbers:
            # Check A + F = 8
            a = 8 - f
            if a in numbers:
                # Now check L - A = 33
                l = a + 33
                if l in numbers:
                    # Check F + L = 41
                    if f + l == 41:
                        # We found F, K, A, L
                        result[0] = a  # A
                        result[10] = k  # K
                        result[5] = int(f)  # F
                        result[11] = l  # L
                        
                        # From E = 2.8H and K - E = -26
                        e = k + 26
                        if e in numbers:
                            h = e / 2.8
                            if h in numbers:
                                if abs(h - 2.0*f) < 0.1:  # H = 2.0F
                                    result[4] = e  # E
                                    result[7] = int(h)  # H
                                    
                                    # From I - E = -4
                                    i = e - 4
                                    if i in numbers:
                                        # Check I = 1.5B
                                        b = i / 1.5
                                        if b in numbers:
                                            result[8] = int(i)  # I
                                            result[1] = int(b)  # B
                                            
                                            # From L - C = -60
                                            c = l + 60
                                            if c in numbers:
                                                # Check C + M = 111
                                                m = 111 - c
                                                if m in numbers:
                                                    result[2] = c  # C
                                                    result[12] = m  # M
                                                    
                                                    # Remaining numbers go to D, G, J
                                                    remaining = numbers - set(result)
                                                    if len(remaining) == 3:
                                                        # D, G, J can be any of remaining numbers
                                                        result[3] = min(remaining)  # D
                                                        remaining.remove(result[3])
                                                        result[6] = min(remaining)  # G
                                                        remaining.remove(result[6])
                                                        result[9] = remaining.pop()  # J
                                                        print(result)
                                                        return

solve_equations()
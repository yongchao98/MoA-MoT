def find_solution():
    # From J = 2.0K, K must be even and J must be double K
    possible_K = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]
    possible_J = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]
    
    for k in possible_K:
        j = 2.0 * k
        if j in possible_J:
            # From L = 1.4K
            l = 1.4 * k
            if l in [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]:
                # From J - G = -86
                g = j + 86
                if g in [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]:
                    remaining = set([2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96])
                    remaining -= set([j, k, l, g])
                    
                    # Try each remaining number for I (must be > L)
                    for i in remaining:
                        if i > l:
                            # From F - I = -13
                            f = i - 13
                            if f in remaining:
                                temp_remaining = remaining - set([i, f])
                                
                                # From L - A = -29
                                a = l + 29
                                if a in temp_remaining:
                                    new_remaining = temp_remaining - set([a])
                                    
                                    # Try each remaining number for B
                                    for b in new_remaining:
                                        # From B + K = 50
                                        if b + k == 50:
                                            # From M - B = -21
                                            m = b - 21
                                            if m in new_remaining:
                                                # From C + M = 52
                                                c = 52 - m
                                                if c in new_remaining:
                                                    # From B + C = 73
                                                    if b + c == 73:
                                                        final_remaining = new_remaining - set([b, m, c])
                                                        
                                                        # Try remaining combinations for D, H
                                                        for d in final_remaining:
                                                            h = 53 - d
                                                            if h in final_remaining:
                                                                # From D + M = 74
                                                                if d + m == 74:
                                                                    e = (final_remaining - set([d, h])).pop()
                                                                    result = [a, b, c, d, e, f, g, h, i, j, k, l, m]
                                                                    print(result)
                                                                    return

find_solution()
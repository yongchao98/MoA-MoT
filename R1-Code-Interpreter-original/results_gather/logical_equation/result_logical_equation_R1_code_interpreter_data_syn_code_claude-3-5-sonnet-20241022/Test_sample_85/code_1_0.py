def find_solution():
    numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]
    
    # From condition 9: C = 1.5K
    # Since K must be an integer from our list, and C must also be from our list
    # Let's find possible K and C pairs first
    possible_K_C = []
    for k in numbers:
        c = 1.5 * k
        if c in numbers:
            possible_K_C.append((k, int(c)))

    # From conditions 4 and 5: D + H = 66 and L - D = -35
    # This means D + H = 66 and L = D - 35
    possible_D_H_L = []
    for d in numbers:
        for h in numbers:
            if d + h == 66:
                l = d - 35
                if l in numbers:
                    possible_D_H_L.append((d, h, l))

    # From condition 2: G = 4.8M
    possible_G_M = []
    for m in numbers:
        g = 4.8 * m
        if abs(g - round(g)) < 0.1 and round(g) in numbers:
            possible_G_M.append((round(g), m))

    # From condition 6: E = 3.6B
    possible_E_B = []
    for b in numbers:
        e = 3.6 * b
        if abs(e - round(e)) < 0.1 and round(e) in numbers:
            possible_E_B.append((round(e), b))

    # Try combinations of our found pairs
    for k, c in possible_K_C:
        for d, h, l in possible_D_H_L:
            for g, m in possible_G_M:
                for e, b in possible_E_B:
                    # Check all conditions
                    if (e > c and  # 1
                        h - g == -8 and  # 8
                        b - d == -40 and  # 7
                        abs(h - 1.6 * b) < 0.1):  # 10
                        
                        # Find J from condition 11: J - E = 9
                        j = e + 9
                        
                        # Check if J is in our numbers
                        if j in numbers:
                            # Find remaining number for A, F, I
                            remaining = set(numbers) - {b, c, d, e, g, h, j, k, l, m}
                            if len(remaining) == 3:
                                # We found our solution
                                result = [0] * 13
                                # A, B, C, D, E, F, G, H, I, J, K, L, M
                                positions = {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 
                                           'F': 5, 'G': 6, 'H': 7, 'I': 8, 'J': 9, 
                                           'K': 10, 'L': 11, 'M': 12}
                                
                                result[positions['B']] = b
                                result[positions['C']] = c
                                result[positions['D']] = d
                                result[positions['E']] = e
                                result[positions['G']] = g
                                result[positions['H']] = h
                                result[positions['J']] = j
                                result[positions['K']] = k
                                result[positions['L']] = l
                                result[positions['M']] = m
                                
                                remaining = list(remaining)
                                result[positions['A']] = remaining[0]
                                result[positions['F']] = remaining[1]
                                result[positions['I']] = remaining[2]
                                
                                print(result)
                                return

find_solution()
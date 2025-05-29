def find_solution():
    # Given numbers
    numbers = {1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50}
    
    # From D = 3.0I, possible pairs of (D,I):
    d_i_pairs = [(d, d/3) for d in numbers if d/3 in numbers]
    
    # From D + H = 8, possible pairs of (D,H):
    d_h_pairs = [(d, h) for d in numbers for h in numbers if d + h == 8]
    
    # From C + H = 20, possible pairs of (C,H):
    c_h_pairs = [(c, h) for c in numbers for h in numbers if c + h == 20]
    
    # From F + G = 38, possible pairs of (F,G):
    f_g_pairs = [(f, g) for f in numbers for g in numbers if f + g == 38]
    
    # From L - G = -26, possible pairs of (L,G):
    l_g_pairs = [(l, g) for l in numbers for g in numbers if l - g == -26]
    
    # From C + F = 17, possible pairs of (C,F):
    c_f_pairs = [(c, f) for c in numbers for f in numbers if c + f == 17]
    
    # From E + L = 60, possible pairs of (E,L):
    e_l_pairs = [(e, l) for e in numbers for l in numbers if e + l == 60]
    
    # From M = 4.5L, possible pairs of (M,L):
    m_l_pairs = [(m, l) for m in numbers for l in numbers if m == 4.5 * l]
    
    # From A - K = -17, possible pairs of (A,K):
    a_k_pairs = [(a, k) for a in numbers for k in numbers if a - k == -17]
    
    for d, i in d_i_pairs:
        for h in numbers:
            if d + h != 8: continue
            for c in numbers:
                if c + h != 20: continue
                for f in numbers:
                    if c + f != 17: continue
                    for g in numbers:
                        if f + g != 38: continue
                        for l in numbers:
                            if l - g != -26: continue
                            for e in numbers:
                                if e + l != 60: continue
                                for m in numbers:
                                    if m != 4.5 * l: continue
                                    for a in numbers:
                                        for k in numbers:
                                            if a - k != -17: continue
                                            if a <= f: continue
                                            
                                            # Check if we used each number exactly once
                                            used = {a,c,d,e,f,g,h,i,k,l,m}
                                            remaining = numbers - used
                                            if len(remaining) == 2:  # For B and J
                                                b = max(remaining)
                                                j = min(remaining)
                                                result = [a,b,c,d,e,f,g,h,i,j,k,l,m]
                                                if len(set(result)) == 13:
                                                    print(result)
                                                    return

find_solution()
def find_solution():
    numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]
    
    # Find possible B values (H = 4.8B must give a value in our list)
    possible_B = []
    for b in numbers:
        h = 4.8 * b
        if abs(h - round(h)) < 0.1 and round(h) in numbers:
            possible_B.append(b)
    
    # For each possible B
    for B in possible_B:
        H = int(4.8 * B)
        
        # From B - I = 3, find I
        I = B - 3
        if I not in numbers:
            continue
            
        # From I = 2F, find F
        F = I / 2
        if F not in numbers:
            continue
            
        # From A - F = 2, find A
        A = F + 2
        if A not in numbers:
            continue
            
        # From A - G = -13, find G
        G = A + 13
        if G not in numbers:
            continue
            
        # From B + C = 15, find C
        C = 15 - B
        if C not in numbers:
            continue
            
        # From C - H = -14, verify
        if C - H != -14:
            continue
            
        # From L - J = -5, try possible L,J combinations
        remaining = set(numbers) - {A, B, C, F, G, H, I}
        for J in remaining:
            L = J - 5
            if L in remaining and L != J:
                # Found a valid solution
                # Last number is M
                M = (set(numbers) - {A, B, C, F, G, H, I, J, L}).pop()
                
                # Create result in alphabetical order
                result = [A, B, C, None, None, F, G, H, I, J, None, L, M]
                
                # Fill in remaining values
                unused = set(numbers) - set(x for x in result if x is not None)
                D, E, K = sorted(unused)
                result[3] = D
                result[4] = E
                result[10] = K
                
                print(result)
                return

find_solution()
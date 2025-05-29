def find_solution():
    numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]
    
    # We can deduce some relationships:
    # 1. From B = 3.2C: if C = 5, then B = 16
    # 2. From H - B = 34: if B = 16, then H = 50
    # 3. From H + L = 60: if H = 50, then L = 10
    # 4. From I = 2.4F: if F = 15, then I = 36
    
    # Start with these known values
    C = 5
    B = 16
    H = 50
    L = 10
    F = 15
    I = 36
    
    # From D = 4.8C: D should be 24
    D = 24
    
    # From G - J = -17: try remaining numbers for G and J
    remaining = set(numbers) - {C, B, H, L, F, I, D}
    
    for G in remaining:
        J = G + 17
        if J not in remaining:
            continue
            
        remaining_2 = remaining - {G, J}
        
        # M = 3.5E
        for E in remaining_2:
            M = int(3.5 * E)
            if M not in remaining_2 or M == E:
                continue
                
            remaining_3 = remaining_2 - {E, M}
            
            # Last number must be A, and H > A must be satisfied
            A = list(remaining_3)[0]
            if H <= A:
                continue
                
            # Check if all conditions are satisfied
            if (I == 2.4 * F and 
                D == 4.8 * C and 
                D == 1.6 * F and 
                B == 1.6 * L and 
                M == 3.5 * E and 
                G - J == -17 and 
                H - B == 34 and 
                B == 3.2 * C and 
                H > E and 
                H + L == 60 and 
                I - C == 31 and 
                H > A):
                
                return [A, B, C, D, E, F, G, H, I, J, K, L, M]
    
    return None

result = find_solution()
if result:
    print(result)
def find_solution():
    numbers = {1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50}
    
    # Try possible B and J combinations that sum to 16
    for B in [10, 9, 8]:  # B must be larger
        J = 16 - B
        if J not in numbers or B not in numbers:
            continue
            
        # Find E based on A + B = 22 and A = 1.4E
        for E in [5, 7, 10, 15]:
            A = int(1.4 * E)
            if A not in numbers or A + B != 22:
                continue
                
            # Check F = 2.0E
            F = 2 * E
            if F not in numbers:
                continue
                
            # Try to find H where A - H = 5
            H = A - 5
            if H not in numbers:
                continue
                
            # Find G where G + J = 25
            G = 25 - J
            if G not in numbers:
                continue
                
            # Find C where C + E = 21
            C = 21 - E
            if C not in numbers:
                continue
                
            # Try D and L where D + L = 48
            remaining = numbers - {A, B, C, E, F, G, H, J}
            for D in remaining:
                L = 48 - D
                if L not in remaining:
                    continue
                    
                # The remaining number must be M
                remaining = numbers - {A, B, C, D, E, F, G, H, J, L}
                if len(remaining) == 2:  # One for I, one for K
                    I = max(remaining)  # Assign largest remaining to I
                    K = min(remaining)  # Assign smallest remaining to K
                    M = (numbers - {A, B, C, D, E, F, G, H, I, J, K, L}).pop()
                    
                    result = [A, B, C, D, E, F, G, H, I, J, K, L, M]
                    print(result)
                    return

find_solution()
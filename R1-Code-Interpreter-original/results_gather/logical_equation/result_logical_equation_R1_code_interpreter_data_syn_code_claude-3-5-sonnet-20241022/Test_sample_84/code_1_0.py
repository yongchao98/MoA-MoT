def find_solution():
    # Available numbers
    numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]
    
    # We know several direct relationships:
    # A + J = 31
    # A + B = 23
    # M - B = 89
    # I - B = 43
    
    for B in numbers:
        # M must be B + 89
        M = B + 89
        if M not in numbers:
            continue
            
        # I must be B + 43
        I = B + 43
        if I not in numbers:
            continue
            
        for G in numbers:
            # Check B = 1.4G
            if abs(B - 1.4 * G) > 0.1:
                continue
                
            # Calculate A using A = 3.2G
            A = int(3.2 * G + 0.1)
            if A not in numbers:
                continue
                
            # Calculate J using J = 3.0G
            J = int(3.0 * G + 0.1)
            if J not in numbers:
                continue
                
            # Verify A + J = 31
            if A + J != 31:
                continue
                
            # Verify A + B = 23
            if A + B != 23:
                continue
                
            # Calculate K using A = 1.6K
            K = int(A / 1.6 + 0.1)
            if K not in numbers:
                continue
                
            # Calculate F using F = 2.4K
            F = int(2.4 * K + 0.1)
            if F not in numbers:
                continue
                
            # Verify F = 4.8G
            if abs(F - 4.8 * G) > 0.1:
                continue
                
            # Calculate L using B = 3.5L
            L = int(B / 3.5 + 0.1)
            if L not in numbers:
                continue
                
            # Find remaining numbers for C, D, E, H
            remaining = set(numbers) - {A, B, F, G, H, I, J, K, L, M}
            if len(remaining) == 3:
                C, D, E = sorted(list(remaining))
                result = [A, B, C, D, E, F, G, H, I, J, K, L, M]
                print(result)
                return

find_solution()
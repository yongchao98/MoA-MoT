# Possible values for the letters
values = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Iterate over possible values for A
for A in values:
    # Check if B = 4.5A is a valid value
    B = 4.5 * A
    if B not in values:
        continue
    
    # Calculate K from A + K = 15
    K = 15 - A
    if K not in values:
        continue
    
    # Calculate M from A - M = 7
    M = A - 7
    if M not in values:
        continue
    
    # Calculate L from M = 3.0L
    L = M / 3.0
    if L not in values:
        continue
    
    # Calculate E from B + E = 52
    E = 52 - B
    if E not in values:
        continue
    
    # Calculate H from H = 1.5A
    H = 1.5 * A
    if H not in values:
        continue
    
    # Check remaining constraints
    for F in values:
        J = 86 - F
        if J not in values:
            continue
        
        for D in values:
            if J <= D:
                continue
            
            for G in values:
                if A <= G:
                    continue
                
                if H <= L:
                    continue
                
                # Assign remaining variables
                for C in values:
                    for I in values:
                        # Check if all values are unique
                        if len(set([A, B, C, D, E, F, G, H, I, J, K, L, M])) == 13:
                            print([A, B, C, D, E, F, G, H, I, J, K, L, M])
                            break
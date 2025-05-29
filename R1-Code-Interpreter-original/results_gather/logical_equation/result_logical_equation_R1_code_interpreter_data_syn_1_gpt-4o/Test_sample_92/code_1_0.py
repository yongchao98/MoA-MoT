# Given numbers
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over possible values for L
for L in numbers:
    A = 3.0 * L
    if A not in numbers:
        continue
    
    M = 39 - A
    if M not in numbers:
        continue
    
    if L + M != 29:
        continue
    
    B = 3.2 * L
    if B not in numbers:
        continue
    
    I = A / 1.5
    if I not in numbers:
        continue
    
    E = 2.4 * A
    if E not in numbers:
        continue
    
    H = A - 12
    if H not in numbers:
        continue
    
    # Iterate over possible values for C and D
    for C in numbers:
        D = 124 - C
        if D not in numbers or C <= D:
            continue
        
        if C - B != 80:
            continue
        
        # Iterate over possible values for F and G
        for F in numbers:
            G = 47 - F
            if G not in numbers:
                continue
            
            # Check E + K = 43
            K = 43 - E
            if K not in numbers:
                continue
            
            # Check if all numbers are unique
            if len(set([A, B, C, D, E, F, G, H, I, K, L, M])) == 12:
                J = list(set(numbers) - set([A, B, C, D, E, F, G, H, I, K, L, M]))[0]
                print([A, B, C, D, E, F, G, H, I, J, K, L, M])
                break
# Given numbers
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over possible values for F
for F in numbers:
    G = 46 - F
    if G not in numbers or G != 3.6 * F:
        continue
    
    B = 1.5 * F
    if B not in numbers:
        continue
    
    M = 31 - B
    if M not in numbers:
        continue
    
    J = M + 80
    if J not in numbers:
        continue
    
    L = 1.5 * M
    if L not in numbers:
        continue
    
    I = L + 21
    if I not in numbers:
        continue
    
    for A in numbers:
        C = 4.0 * A
        if C not in numbers:
            continue
        
        for E in numbers:
            D = 1.5 * E
            if D not in numbers:
                continue
            
            H = L + 26
            if H not in numbers:
                continue
            
            # Check if all numbers are used
            if set([A, B, C, D, E, F, G, H, I, J, L, M]) == set(numbers):
                print([A, B, C, D, E, F, G, H, I, J, K, L, M])
                break
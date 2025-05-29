# Given numbers
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over possible values for C
for C in numbers:
    E = C + 11
    if E not in numbers:
        continue
    
    J = 4.8 * C
    if J not in numbers:
        continue
    
    D = J - 14
    if D not in numbers:
        continue
    
    H = J - 22
    if H not in numbers:
        continue
    
    for L in numbers:
        I = L - 12
        M = 2.4 * L
        if I not in numbers or M not in numbers:
            continue
        
        K = 3.5 * H
        if K not in numbers:
            continue
        
        # Check if all values are distinct
        if len(set([C, D, E, H, I, J, K, L, M])) == 9:
            # Find remaining numbers for A, B, F, G
            remaining = set(numbers) - {C, D, E, H, I, J, K, L, M}
            if len(remaining) == 4:
                A, B, F, G = remaining
                # Print the result
                print([A, B, C, D, E, F, G, H, I, J, K, L, M])
                break
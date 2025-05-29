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
    
    I = A / 1.5
    if I not in numbers:
        continue
    
    E = 2.4 * A
    if E not in numbers:
        continue
    
    B = 3.2 * L
    if B not in numbers:
        continue
    
    C = B + 80
    if C not in numbers:
        continue
    
    D = 124 - C
    if D not in numbers or C <= D:
        continue
    
    H = A - 12
    if H not in numbers:
        continue
    
    # Check F + G = 47
    remaining_numbers = set(numbers) - {A, B, C, D, E, H, I, L, M}
    for F in remaining_numbers:
        G = 47 - F
        if G in remaining_numbers and F != G:
            K = 43 - E
            if K in remaining_numbers - {F, G}:
                J = (remaining_numbers - {F, G, K}).pop()
                print([A, B, C, D, E, F, G, H, I, J, K, L, M])
                break
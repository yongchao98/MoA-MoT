# Given numbers
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over possible values for J
for J in numbers:
    L = 3 * J
    if L not in numbers:
        continue
    
    D = L + 30
    if D not in numbers:
        continue
    
    F = 17 - L
    if F not in numbers:
        continue
    
    B = F + 22
    if B not in numbers:
        continue
    
    A = 4 * B
    if A not in numbers:
        continue
    
    I = D - 17
    if I not in numbers:
        continue
    
    E = I - 21
    if E not in numbers:
        continue
    
    C = D - 29
    if C not in numbers:
        continue
    
    # Check if all values are unique
    if len({A, B, C, D, E, F, I, J, L}) == 9:
        # Find remaining numbers for G, H, K, M
        remaining = set(numbers) - {A, B, C, D, E, F, I, J, L}
        G, H, K, M = remaining
        # Print the result
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
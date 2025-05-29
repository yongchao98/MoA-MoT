# Possible values for the letters
values = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Iterate over possible values for J, since A and G depend on it
for J in values:
    A = J + 6
    if A not in values:
        continue
    
    B = A + 9
    if B not in values:
        continue
    
    H = 12 - B
    if H not in values:
        continue
    
    L = B + 17
    if L not in values:
        continue
    
    C = 52 - B
    if C not in values:
        continue
    
    K = C + 5
    if K not in values:
        continue
    
    E = C / 3
    if E not in values:
        continue
    
    D = E - 13
    if D not in values:
        continue
    
    I = 18 - E
    if I not in values:
        continue
    
    G = 3.6 * J
    if G not in values:
        continue
    
    # Check if all values are unique
    if len(set([A, B, C, D, E, F, G, H, I, J, K, L, M])) == 13:
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break
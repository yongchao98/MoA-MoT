def find_solution():
    numbers = set([2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96])
    
    # From I + M = 7, since both are positive integers from our set
    # M must be 2 and I must be 5
    M = 2
    I = 5
    numbers.remove(M)
    numbers.remove(I)
    
    # From F = 3.5M and F = 1.4I
    # F = 7 (3.5 * 2 = 7 and 1.4 * 5 = 7)
    F = 7
    numbers.remove(F)
    
    # From M - H = -48
    # H = M + 48 = 2 + 48 = 50
    H = 50
    numbers.remove(H)
    
    # From A - H = -14
    # A = H - 14 = 50 - 14 = 36
    A = 36
    numbers.remove(A)
    
    # From A - G = 33
    # G = A - 33 = 36 - 33 = 3
    G = 3
    numbers.remove(G)
    
    # From K - G = 25
    # K = G + 25 = 3 + 25 = 28
    K = 28
    numbers.remove(K)
    
    # From J + K = 38
    # J = 38 - K = 38 - 28 = 10
    J = 10
    numbers.remove(J)
    
    # From B + I = 50
    # B = 50 - I = 50 - 5 = 45
    B = 45
    numbers.remove(B)
    
    # From B + E = 141
    # E = 141 - B = 141 - 45 = 96
    E = 96
    numbers.remove(E)
    
    # From C - E = -72
    # C = E - 72 = 96 - 72 = 24
    C = 24
    numbers.remove(C)
    
    # Only remaining numbers are for D and L
    # D = 15 and L = 16 (can be verified)
    D = 15
    L = 16
    
    result = [A, B, C, D, E, F, G, H, I, J, K, L, M]
    print(result)

find_solution()
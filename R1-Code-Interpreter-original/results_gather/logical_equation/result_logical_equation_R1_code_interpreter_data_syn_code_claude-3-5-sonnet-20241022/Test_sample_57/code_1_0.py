def find_solution():
    # Known values from direct equations:
    # C - G = -36 and C + G = 42 implies:
    # 2C = 6, so C = 3
    # G = 39
    
    numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]
    C = 3
    G = 39
    
    # From C - F = -17, we know F = 20
    F = 20
    
    # From B + F = 68, we know B = 48
    B = 48
    
    # From C + H = 124, we know H = 121
    H = 121
    
    # From E + H = 137, we know E = 16
    E = 16
    
    # From A + K = 375 and A > C, A > F
    # A must be 150 and K must be 225
    A = 150
    K = 225
    
    # J > C means J must be > 3
    # Among remaining numbers [9, 75, 80]
    # J = 80 (needs to be greater than C but not too large)
    J = 80
    
    # Remaining numbers: [9, 75]
    # D = 75 (larger remaining number)
    # I = 9 (last remaining number)
    D = 75
    I = 9
    
    result = [A, B, C, D, E, F, G, H, I, J, K]
    
    # Verify all constraints
    assert C - G == -36
    assert C - F == -17
    assert B + F == 68
    assert E + H == 137
    assert A + K == 375
    assert C + H == 124
    assert A > F
    assert J > C
    assert F - G == -19
    assert A > C
    assert C + G == 42
    
    print(result)

find_solution()
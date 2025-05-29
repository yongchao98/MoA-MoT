def find_solution():
    # Given numbers
    nums = {1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50}
    
    # From I = 3.6L and I = 2.4E, we can deduce:
    # I must be 36 (as it's the only number divisible by both 2.4 and 3.6)
    I = 36
    # Therefore L = 10 (36/3.6)
    L = 10
    # And E = 15 (36/2.4)
    E = 15
    
    # From J = 1.6E
    # J = 1.6 * 15 = 24
    J = 24
    
    # From E + H = 31
    # H = 31 - 15 = 16
    H = 16
    
    # From A - H = -14
    # A = H - 14 = 16 - 14 = 2
    A = 2
    
    # From C - I = -35
    # C = I - 35 = 36 - 35 = 1
    C = 1
    
    # From D + K = 31 and F - K = 4
    remaining = nums - {A, C, E, H, I, J, L}
    
    for d in remaining:
        k = 31 - d  # from D + K = 31
        if k not in remaining:
            continue
        f = k + 4   # from F - K = 4
        if f not in remaining:
            continue
            
        # From D + G = 33
        g = 33 - d
        if g not in remaining:
            continue
            
        if g >= J:  # J > G must be satisfied
            continue
            
        # Check if we have valid values
        used = {A, C, d, E, f, g, H, I, J, k, L}
        remaining_nums = nums - used
        if len(remaining_nums) != 2:  # Should have exactly 2 numbers left for B and M
            continue
            
        # The larger number should be M
        m = max(remaining_nums)
        b = min(remaining_nums)
        
        # Return the result in order [A, B, C, D, E, F, G, H, I, J, K, L, M]
        return [A, b, C, d, E, f, g, H, I, J, k, L, m]

result = find_solution()
print(result)
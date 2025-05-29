def verify_all_constraints():
    # Fixed values based on equations
    B = 121  # From H - B = -117 and H = 4
    C = 1    # From C + H = 5 and H = 4C
    D = 9    # From D + H = 13 and H = 4
    E = 16   # From E - C = 15 and C = 1
    F = 3    # From F = 3C and C = 1
    H = 4    # From multiple equations
    
    # Remaining numbers: 27, 36, 80 for positions A, G, I
    # Let's try all possibilities in ascending order
    remaining = [27, 36, 80]
    
    # Try the first valid arrangement
    A = 27
    G = 36
    I = 80
    
    solution = [A, B, C, D, E, F, G, H, I]
    
    # Verify one final time
    print(f"Final verification of solution: {solution}")
    print(f"All numbers used exactly once: {sorted(solution) == sorted([1, 3, 4, 9, 16, 27, 36, 80, 121])}")
    print(f"All equations satisfied:")
    print(f"D + H = 13: {D + H == 13}")
    print(f"E - C = 15: {E - C == 15}")
    print(f"C + H = 5: {C + H == 5}")
    print(f"H = 4.0C: {H == 4.0 * C}")
    print(f"H - B = -117: {H - B == -117}")
    print(f"F = 3.0C: {F == 3.0 * C}")
    print(f"F + H = 7: {F + H == 7}")
    
    return solution

verify_all_constraints()
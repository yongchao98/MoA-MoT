def is_valid_partial(F, H, B, D, E, M, I, J, L):
    # Quick checks that can be done with individual or pairs of numbers
    if F and H and F + H != 6:  # Condition 1
        return False
    if F and B and abs(B - 2.0*F) > 0.01:  # Condition 3
        return False
    if D and E and abs(E - 2.4*D) > 0.01:  # Condition 2
        return False
    if B and I and I - B != 34:  # Condition 5
        return False
    if B and H and abs(H - 2.5*B) > 0.01:  # Condition 6
        return False
    if E and M and E + M != 39:  # Condition 7
        return False
    if J and L and abs(L - 4.0*J) > 0.01:  # Condition 8
        return False
    if H and M and abs(M - 3.0*H) > 0.01:  # Condition 9
        return False
    if D and L and abs(L - 2.8*D) > 0.01:  # Condition 10
        return False
    return True

# Start with F and H since they must sum to 6
numbers = {1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50}
result = [0] * 13

for F in [1, 2, 3, 5]:  # F must be small enough for F + H = 6
    H = 6 - F
    if H not in numbers:
        continue
        
    # B must be 2*F
    B = 2 * F
    if B not in numbers:
        continue
        
    # I must be B + 34
    I = B + 34
    if I not in numbers:
        continue
        
    # H must be 2.5*B
    if abs(H - 2.5*B) > 0.01:
        continue
        
    # M must be 3*H
    M = 3 * H
    if M not in numbers:
        continue
        
    remaining = numbers - {F, H, B, I, M}
    
    # Try values for D
    for D in remaining:
        E = int(2.4 * D)
        if E not in remaining or E + M != 39:
            continue
            
        rem2 = remaining - {D, E}
        
        # Try values for J
        for J in rem2:
            L = 4 * J
            if L not in rem2 or abs(L - 2.8*D) > 0.01:
                continue
                
            rem3 = rem2 - {J, L}
            
            # C must be 3*M
            C = 3 * M
            if C not in rem3:
                continue
                
            # Last remaining number is G
            rem4 = rem3 - {C}
            if len(rem4) == 1:
                G = rem4.pop()
                # We found the solution
                result = [7, B, C, D, E, F, G, H, I, J, 24, L, M]
                print(result)
                exit()
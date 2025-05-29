def solve():
    # Possible values for the letters
    values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]
    
    # Precompute possible values for K, J, and L based on the equations
    possible_K = [v for v in values if 2.0 * v in values and 1.4 * v in values]
    
    for K in possible_K:
        J = 2.0 * K
        L = 1.4 * K
        
        if J not in values or L not in values:
            continue
        
        # Remove K, J, L from the list of available values
        remaining_values = [v for v in values if v not in {K, J, L}]
        
        # Use constraints to find B and M
        for B in remaining_values:
            if 50 - B in remaining_values:
                K = 50 - B
                if K not in remaining_values:
                    continue
                
                M = B + 21
                if M not in remaining_values:
                    continue
                
                # Remove B, K, M from the list of available values
                remaining_values_2 = [v for v in remaining_values if v not in {B, K, M}]
                
                # Use remaining constraints to find other values
                for C in remaining_values_2:
                    if 73 - C in remaining_values_2:
                        D = 74 - M
                        if D not in remaining_values_2:
                            continue
                        
                        H = 53 - D
                        if H not in remaining_values_2:
                            continue
                        
                        G = J - 86
                        if G not in remaining_values_2:
                            continue
                        
                        for F in remaining_values_2:
                            if F - 13 in remaining_values_2:
                                I = F - 13
                                if I not in remaining_values_2:
                                    continue
                                
                                if I <= L:
                                    continue
                                
                                # Find remaining values for A, E
                                remaining_values_3 = [v for v in remaining_values_2 if v not in {C, D, H, G, F, I}]
                                
                                for A in remaining_values_3:
                                    if L - A == 29:
                                        E = [v for v in remaining_values_3 if v != A][0]
                                        
                                        # Print the result in alphabetical order
                                        print(f"<<<[{A}, {B}, {C}, {D}, {E}, {F}, {G}, {H}, {I}, {J}, {K}, {L}, {M}]>>>")
                                        return

solve()
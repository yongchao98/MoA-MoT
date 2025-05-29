from itertools import permutations

# Possible values for the letters
values = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Iterate over all possible values for E
for E in values:
    F = 2.0 * E
    A = 1.4 * E
    
    # Check if F and A are in the list of possible values
    if F not in values or A not in values:
        continue
    
    # Iterate over all possible values for B
    for B in values:
        if B <= E:
            continue
        
        J = 16 - B
        if J not in values or B <= J:
            continue
        
        G = 25 - J
        if G not in values:
            continue
        
        # Iterate over all possible values for D
        for D in values:
            L = 48 - D
            if L not in values:
                continue
            
            # Iterate over all possible values for C
            for C in values:
                if C + E != 21:
                    continue
                
                # Calculate H from A - H = 5
                H = A - 5
                if H not in values:
                    continue
                
                # Check if all values are unique
                if len(set([A, B, C, D, E, F, G, H, J, L])) == 10:
                    # Find remaining values for I, K, M
                    remaining_values = set(values) - set([A, B, C, D, E, F, G, H, J, L])
                    I, K, M = remaining_values
                    
                    # Print the solution
                    print([A, B, C, D, E, F, G, H, I, J, K, L, M])
                    break
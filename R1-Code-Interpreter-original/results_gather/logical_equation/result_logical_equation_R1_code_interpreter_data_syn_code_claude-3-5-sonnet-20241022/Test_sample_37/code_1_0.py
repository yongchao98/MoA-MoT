def find_solution():
    numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]
    
    # Pre-filter based on known relationships
    for H in numbers:
        I = 1.5 * H
        if I not in numbers:
            continue
            
        for B in numbers:
            F = 4.0 * B
            if F not in numbers:
                continue
                
            for J in numbers:
                if J - F != -20:  # Constraint 9
                    continue
                    
                G = J + 11  # From constraint 5
                if G not in numbers:
                    continue
                    
                E = 3.0 * J  # From constraint 4
                if E not in numbers:
                    continue
                    
                if F + H != 186:  # Constraint 10
                    continue
                    
                K = 305 - I  # From constraint 7
                if K not in numbers:
                    continue
                    
                D = I - 222  # From constraint 8
                if D not in numbers:
                    continue
                    
                remaining_nums = set(numbers) - {B, D, E, F, G, H, I, J, K}
                for C in remaining_nums:
                    if C + E != 169:  # Constraint 3
                        continue
                    if C <= F:  # Constraint 6
                        continue
                        
                    remaining = list(remaining_nums - {C})
                    if len(remaining) == 1:
                        A = remaining[0]
                        result = [A, B, C, D, E, F, G, H, I, J, K]
                        # Verify all numbers are used exactly once
                        if sorted(result) == sorted(numbers):
                            print(result)
                            return

find_solution()
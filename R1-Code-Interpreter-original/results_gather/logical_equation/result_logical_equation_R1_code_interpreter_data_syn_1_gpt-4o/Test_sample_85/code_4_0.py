from itertools import permutations

# List of numbers to assign
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over possible values for D, K, M
for D in numbers:
    B = D + 40
    L = D - 35
    if B not in numbers or L not in numbers:
        continue
    
    for K in numbers:
        C = 1.5 * K
        if C not in numbers:
            continue
        
        for M in numbers:
            G = 4.8 * M
            if G not in numbers:
                continue
            
            # Calculate other dependent variables
            H = 1.6 * B
            E = 3.6 * B
            if H not in numbers or E not in numbers:
                continue
            
            # Check remaining conditions
            for perm in permutations(numbers):
                A, F, I, J = perm[:4]  # Select 4 numbers for A, F, I, J
                
                if (E > C and
                    L - C == 12 and
                    D + H == 66 and
                    J - E == 9 and
                    H - G == -8):
                    
                    # If all conditions are satisfied, print the result
                    print([A, B, C, D, E, F, G, H, I, J, K, L, M])
                    break
# Given values
values = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Assign known values based on reasoning
A = 10
B = 45
K = 5
M = 3
L = 1
E = 7
H = 15

# Remaining values to assign: C, D, F, G, I, J
remaining_values = set(values) - {A, B, K, M, L, E, H}

# Iterate over possible values for F and J
for F in remaining_values:
    J = 86 - F
    if J not in remaining_values:
        continue
    
    # Iterate over possible values for D, G, I
    for D in remaining_values - {F, J}:
        if J <= D:
            continue
        
        for G in remaining_values - {F, J, D}:
            if A <= G:
                continue
            
            for I in remaining_values - {F, J, D, G}:
                # Assign C as the remaining value
                C = (remaining_values - {F, J, D, G, I}).pop()
                # Check if all values are unique
                if len(set([A, B, C, D, E, F, G, H, I, J, K, L, M])) == 13:
                    print([A, B, C, D, E, F, G, H, I, J, K, L, M])
                    break
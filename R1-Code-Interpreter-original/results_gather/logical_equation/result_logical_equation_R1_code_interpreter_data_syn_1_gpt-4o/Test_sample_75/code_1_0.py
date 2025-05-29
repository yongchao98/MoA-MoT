# Possible values for the letters
values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Find possible values for F and J based on F - J = 38
possible_FJ = [(F, J) for F in values for J in values if F - J == 38]

# Iterate over possible F and J values
for F, J in possible_FJ:
    # Calculate H based on H - F = 5
    H = F + 5
    if H not in values:
        continue
    
    # Calculate C based on C + F = 55
    C = 55 - F
    if C not in values:
        continue
    
    # Calculate I based on F + I = 73
    I = 73 - F
    if I not in values:
        continue
    
    # Calculate M based on I + M = 44
    M = 44 - I
    if M not in values:
        continue
    
    # Calculate L based on L - I = -23
    L = I - 23
    if L not in values:
        continue
    
    # Calculate B based on B - C = 26
    B = C + 26
    if B not in values:
        continue
    
    # Check D + K = 98
    for D in values:
        K = 98 - D
        if K in values and K != D:
            # Calculate A based on B = 2.4A
            if B % 2.4 == 0:
                A = B / 2.4
                if A in values:
                    # Check if all values are unique
                    used_values = {A, B, C, D, E, F, G, H, I, J, K, L, M}
                    if len(used_values) == 13:
                        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
                        break
from sympy import symbols, Eq, solve

# Define the variables
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the possible values
possible_values = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Solve some equations directly
# 1. I = 1.5H
# 2. F = 4.0B
# 3. E = 3.0J
# 4. J - G = -11
# 5. I + K = 305
# 6. I - D = 222
# 7. J - F = -20
# 8. F + H = 186

# Start solving incrementally
for H in possible_values:
    I = 1.5 * H
    if I not in possible_values:
        continue
    
    for B in possible_values:
        F = 4.0 * B
        if F not in possible_values:
            continue
        
        for J in possible_values:
            E = 3.0 * J
            if E not in possible_values:
                continue
            
            G = J + 11
            if G not in possible_values:
                continue
            
            K = 305 - I
            if K not in possible_values:
                continue
            
            D = I - 222
            if D not in possible_values:
                continue
            
            if J - F != -20:
                continue
            
            if F + H != 186:
                continue
            
            C = 169 - E
            if C not in possible_values:
                continue
            
            if C <= F:
                continue
            
            # Check if all values are unique
            used_values = {H, I, B, F, J, E, G, K, D, C}
            if len(used_values) == 10:
                A = list(set(possible_values) - used_values)[0]
                result = [A, B, C, D, E, F, G, H, I, J, K]
                print(f"<<<{result}>>>")
                break
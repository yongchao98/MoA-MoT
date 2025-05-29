from sympy import symbols, Eq, solve

# Define the symbolic variables
A, B, C, D, E, F, G, H, I = symbols('A B C D E F G H I', integer=True)

# Possible values for the letters
possible_values = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Step-by-step solving
# 1. Solve A = C - 44
# 2. Solve B = C - 79
# 3. Solve D = 3 * B
# 4. Solve I = 3 * D
# 5. Solve A = 4 * I
# 6. Solve B = E + 3
# 7. C > H

# Iterate over possible values for C
for C_val in possible_values:
    A_val = C_val - 44
    B_val = C_val - 79
    D_val = 3 * B_val
    I_val = 3 * D_val
    E_val = B_val - 3
    
    # Check if A, B, D, I, E are in possible values
    if (A_val in possible_values and
        B_val in possible_values and
        D_val in possible_values and
        I_val in possible_values and
        E_val in possible_values):
        
        # Check C > H
        for H_val in possible_values:
            if C_val > H_val:
                # Assign remaining values to F, G, H
                remaining_values = set(possible_values) - {A_val, B_val, C_val, D_val, E_val, I_val}
                if len(remaining_values) == 3:
                    F_val, G_val, H_val = remaining_values
                    result = [A_val, B_val, C_val, D_val, E_val, F_val, G_val, H_val, I_val]
                    print(f"<<<{result}>>>")
                    break
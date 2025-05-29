# Possible values for A, B, C, D, E, F, G, H, I, J, K, L, M
possible_values = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Manually solve the equations
# 1. A + H = 25 and H = 1.5A => A + 1.5A = 25 => 2.5A = 25 => A = 10, H = 15
A_val = 10
H_val = 15

# 2. E = 2.4A => E = 2.4 * 10 = 24
E_val = 24

# 3. L = 1.5E => L = 1.5 * 24 = 36
L_val = 36

# 4. K = 2.5C => K must be a multiple of 2.5
# Possible values for C: 1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50
# Check which value of C gives a valid K
for C_val in possible_values:
    if 2.5 * C_val in possible_values:
        K_val = 2.5 * C_val
        break

# 5. D = K + 23
D_val = K_val + 23

# 6. I = G - 42 and I = E - 21 => G = I + 42 and I = 24 - 21 = 3
I_val = 3
G_val = I_val + 42

# Remaining numbers to assign: 1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50
# Assigned numbers: A_val, E_val, H_val, L_val, K_val, C_val, D_val, I_val, G_val

# Find the remaining numbers
assigned_numbers = {A_val, E_val, H_val, L_val, K_val, C_val, D_val, I_val, G_val}
remaining_numbers = list(set(possible_values) - assigned_numbers)

# Assign remaining numbers to B, F, J, M ensuring B > F
B_val, F_val, J_val, M_val = sorted(remaining_numbers, reverse=True)

# Ensure B > F
if B_val < F_val:
    B_val, F_val = F_val, B_val

# Return the values in alphabetical order
result = [A_val, B_val, C_val, D_val, E_val, F_val, G_val, H_val, I_val, J_val, K_val, L_val, M_val]
print(result)
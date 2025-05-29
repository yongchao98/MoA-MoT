from sympy import symbols, Eq, solve

# Define the symbolic variables
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations based on the constraints
equations = [
    Eq(C + H, 52),
    Eq(G + I, 13),
    Eq(H, 2.4 * E),
    Eq(F + J, 57),
    Eq(A - I, 25),
    Eq(J - I, 4),
    Eq(D + G, 34),
    Eq(F + L, 95),
    Eq(A, 4.0 * J),
    Eq(F + G, 60),
    Eq(J, 1.4 * K)
]

# Solve the system of equations for a subset of variables
partial_solution = solve(equations[:6], (A, I, J, K, E, H), dict=True)

# Possible values for the letters
possible_values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Check which partial solution can be extended to a full solution
for sol in partial_solution:
    # Assign known values
    A_val = sol[A]
    I_val = sol[I]
    J_val = sol[J]
    K_val = sol[K]
    E_val = sol[E]
    H_val = sol[H]
    
    # Deduce other values
    G_val = 13 - I_val
    D_val = 34 - G_val
    F_val = 60 - G_val
    L_val = 95 - F_val
    C_val = 52 - H_val
    
    # Check if all deduced values are in the possible values
    values = [A_val, None, C_val, D_val, E_val, F_val, G_val, H_val, I_val, J_val, K_val, L_val, None]
    if all(v in possible_values for v in values if v is not None):
        # Find remaining values for B and M
        remaining_values = set(possible_values) - set(values)
        if len(remaining_values) == 2:
            B_val, M_val = remaining_values
            values[1] = B_val
            values[12] = M_val
            print(f"<<<{values}>>>")
            break
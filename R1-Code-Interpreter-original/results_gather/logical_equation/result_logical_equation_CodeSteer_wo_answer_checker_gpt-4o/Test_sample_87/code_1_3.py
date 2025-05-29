from sympy import symbols, Eq, solve

# Define symbolic variables for each letter
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations based on the given inequalities
equations = [
    Eq(I - J, -21),
    Eq(J - D, 14),
    Eq(K, 3.5 * H),
    Eq(E, 3.2 * C),
    Eq(J, 4.8 * C),
    Eq(E, 1.6 * D),
    Eq(I - L, -12),
    Eq(M, 2.4 * L),
    Eq(C - E, -11)
]

# Solve the system of equations
solutions = solve(equations, (C, D, E, H, I, J, K, L, M), dict=True)

# Given numbers
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Extract the solutions and map them to the given numbers
# We need to ensure that each letter gets a unique value from the numbers
if solutions:
    solution = solutions[0]
    C_val = solution[C]
    D_val = solution[D]
    E_val = solution[E]
    H_val = solution[H]
    I_val = solution[I]
    J_val = solution[J]
    K_val = solution[K]
    L_val = solution[L]
    M_val = solution[M]

    # Used numbers
    used_numbers = {C_val, D_val, E_val, H_val, I_val, J_val, K_val, L_val, M_val}
    remaining_numbers = sorted(set(numbers) - used_numbers)

    # Assign remaining numbers to A, B, F, G
    A_val, B_val, F_val, G_val = remaining_numbers

    # Output the result in alphabetical order
    result = [A_val, B_val, C_val, D_val, E_val, F_val, G_val, H_val, I_val, J_val, K_val, L_val, M_val]
    print(f"<<<{result}>>>")
else:
    print("No solution found.")
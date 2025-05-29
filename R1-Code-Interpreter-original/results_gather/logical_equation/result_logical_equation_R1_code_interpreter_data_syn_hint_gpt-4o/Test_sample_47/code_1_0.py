from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, E, F, H, I, J, K = symbols('A B C E F H I J K')

# Define the equations based on the constraints
eq1 = Eq(I, 4 * C)
eq2 = Eq(F + I, 83)
eq3 = Eq(K - E, 82)
eq4 = Eq(E + F, 42)
eq5 = Eq(A + H, 225)
eq6 = Eq(B + K, 130)
eq7 = Eq(J, 3 * H)
eq8 = Eq(B, 3 * F)
eq9 = Eq(C + F, 23)

# Solve the system of equations
solution = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9), (A, B, C, E, F, H, I, J, K), dict=True)

# Check if a solution exists
if solution:
    # Extract the solution
    solution = solution[0]
    A_val = solution[A]
    B_val = solution[B]
    C_val = solution[C]
    E_val = solution[E]
    F_val = solution[F]
    H_val = solution[H]
    I_val = solution[I]
    J_val = solution[J]
    K_val = solution[K]

    # Print the values in alphabetical order
    print([A_val, B_val, C_val, E_val, F_val, H_val, I_val, J_val, K_val])
else:
    print("No solution found")
from sympy import symbols, Eq, nonlinsolve, S

# Define the symbolic variables
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M', integer=True)

# Define the equations
equations = [
    Eq(C + E, 141),
    Eq(J + L, 33),
    Eq(K, 2.4 * H),
    Eq(B - J, 10),
    Eq(E - C, 51),
    Eq(E + K, 120),
    Eq(J - H, -5),
    Eq(K - B, 9),
    Eq(B - M, -1),
    Eq(M, 3.2 * J)
]

# Define the possible values for each letter
possible_values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Solve the system of equations
solutions = nonlinsolve(equations, (A, B, C, D, E, F, G, H, I, J, K, L, M))

# Filter the solutions
for sol in solutions:
    # Convert the solution to a dictionary
    sol_dict = {var: sol[i].evalf() for i, var in enumerate((A, B, C, D, E, F, G, H, I, J, K, L, M))}
    
    # Check if the solution satisfies the inequalities and possible values
    if (sol_dict[C] > sol_dict[F] and sol_dict[C] > sol_dict[M] and
        all(int(sol_dict[var]) == sol_dict[var] and sol_dict[var] in possible_values for var in sol_dict)):
        
        # Create the result list in alphabetical order
        result = [int(sol_dict[A]), int(sol_dict[B]), int(sol_dict[C]), int(sol_dict[D]), int(sol_dict[E]), int(sol_dict[F]),
                  int(sol_dict[G]), int(sol_dict[H]), int(sol_dict[I]), int(sol_dict[J]), int(sol_dict[K]), int(sol_dict[L]), int(sol_dict[M])]
        break

# Print the result in the required format
print(f"<<<{result}>>>")
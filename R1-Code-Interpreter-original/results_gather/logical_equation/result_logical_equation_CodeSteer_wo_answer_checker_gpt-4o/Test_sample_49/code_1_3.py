from sympy import symbols, Eq, nonlinsolve

# Define the symbolic variables
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations based on the problem statement
equations = [
    Eq(A + I, 300),
    Eq(E, 2.4 * C),
    Eq(G, 3.0 * K),
    Eq(A - G, 66),
    Eq(J, 4.0 * C),
    Eq(I, 3.0 * A),
    Eq(F, 2.0 * A),
    Eq(I, 1.5 * F),
    Eq(A - C, 55),
    Eq(H - A, -59)
]

# Define the possible values for each letter
possible_values = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Solve the system of equations
solutions = nonlinsolve(equations, [A, C, E, F, G, H, I, J, K])

# Check each solution against the inequalities and possible values
for sol in solutions:
    sol_dict = {A: sol[0], C: sol[1], E: sol[2], F: sol[3], G: sol[4], H: sol[5], I: sol[6], J: sol[7], K: sol[8]}
    if sol_dict[J] > sol_dict[A] and all(value in possible_values for value in sol_dict.values()):
        # Assign B and D to the remaining values
        remaining_values = set(possible_values) - set(sol_dict.values())
        if len(remaining_values) == 2:
            sol_dict[B], sol_dict[D] = remaining_values
            result = [sol_dict[A], sol_dict[B], sol_dict[C], sol_dict[D], sol_dict[E], sol_dict[F], sol_dict[G], sol_dict[H], sol_dict[I], sol_dict[J], sol_dict[K]]
            print(f"<<<{result}>>>")
            break
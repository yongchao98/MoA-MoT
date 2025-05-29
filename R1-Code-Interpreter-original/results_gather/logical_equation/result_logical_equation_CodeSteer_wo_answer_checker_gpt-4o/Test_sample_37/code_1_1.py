from sympy import symbols, Eq, solve

# Define the variables
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations and inequalities
equations = [
    Eq(I, 1.5 * H),
    Eq(F, 4.0 * B),
    Eq(C + E, 169),
    Eq(E, 3.0 * J),
    Eq(J - G, -11),
    Eq(I + K, 305),
    Eq(I - D, 222),
    Eq(J - F, -20),
    Eq(F + H, 186),
    C > F
]

# Define the possible values
possible_values = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Solve the system
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K), dict=True)

# Filter solutions that use all the possible values
for sol in solution:
    if all(value in possible_values for value in sol.values()):
        result = [sol[A], sol[B], sol[C], sol[D], sol[E], sol[F], sol[G], sol[H], sol[I], sol[J], sol[K]]
        print(f"<<<{result}>>>")
        break
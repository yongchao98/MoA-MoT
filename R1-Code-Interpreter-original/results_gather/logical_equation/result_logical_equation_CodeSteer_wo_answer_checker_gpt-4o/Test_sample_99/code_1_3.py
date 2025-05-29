from sympy import symbols, Eq, nonlinsolve

# Define the symbolic variables
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M', integer=True)

# Define the equations
equations = [
    Eq(B + I, 31),
    Eq(I, 1.5 * D),
    Eq(K - E, -2),
    Eq(B + H, 57),
    Eq(C + F, 132),
    Eq(K - D, -13),
    Eq(M, 2.0 * E),
    Eq(A + G, 43),
    Eq(D - H, -34),
    Eq(F, 4.0 * I),
    Eq(I + L, 69)
]

# Define the possible values
possible_values = {2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96}

# Solve the system of equations
solutions = nonlinsolve(equations, (A, B, C, D, E, F, G, H, I, J, K, L, M))

# Debug: Print the raw solutions
print("Raw solutions:", solutions)

# Filter solutions to match possible values and satisfy the inequality M > B
for sol in solutions:
    # Convert the solution to a list of integers
    sol_values = [int(val) for val in sol]
    
    # Check if all values are in the possible values set and if M > B
    if all(val in possible_values for val in sol_values) and sol[M] > sol[B]:
        result = [sol[A], sol[B], sol[C], sol[D], sol[E], sol[F], sol[G], sol[H], sol[I], sol[J], sol[K], sol[L], sol[M]]
        print(f"<<<{result}>>>")
        break
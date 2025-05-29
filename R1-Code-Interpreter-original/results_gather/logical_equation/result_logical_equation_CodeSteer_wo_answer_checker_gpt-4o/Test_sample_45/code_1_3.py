from sympy import symbols, Eq, nonlinsolve

# Define symbolic variables
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K', integer=True)

# Define the equations and inequalities
equations = [
    Eq(F, 2.4 * D),
    Eq(G, 1.5 * K),
    Eq(G + K, 375),
    Eq(F - G, -177),
    Eq(A + D, 29),
    Eq(I - D, 19),
    Eq(A, 3.0 * H),
    Eq(F, 3.0 * B),
    Eq(F + K, 198),
    Eq(J - F, 73),
    G > C
]

# Given numbers
numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Solve the system of equations and inequalities
solutions = nonlinsolve(equations, (A, B, C, D, E, F, G, H, I, J, K))

# Filter solutions that match the given numbers
for sol in solutions:
    values = list(sol)
    if sorted(values) == sorted(numbers):
        result = values
        break

print(f"<<<{result}>>>")
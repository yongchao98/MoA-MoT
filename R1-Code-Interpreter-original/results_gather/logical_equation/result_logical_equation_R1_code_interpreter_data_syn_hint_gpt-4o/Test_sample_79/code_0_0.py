from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations
eq1 = Eq(F + J, 15)
eq2 = Eq(F - G, -91)
eq3 = Eq(D, 3.5 * K)
eq4 = Eq(E + G, 141)
eq5 = Eq(A - K, 13)
eq6 = Eq(F + M, 29)
eq7 = Eq(D, 1.4 * F)
eq8 = Eq(L, 4.0 * D)
eq9 = Eq(C, 3.6 * J)
eq10 = Eq(A + M, 39)
eq11 = Eq(C, 2.4 * A)

# Solve the system of equations
solution = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11), (A, B, C, D, E, F, G, H, I, J, K, L, M))

# Print the solution
print(solution)
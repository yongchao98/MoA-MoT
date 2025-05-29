from sympy import symbols, Eq, solve

# Define the symbols
L, M, F, H, K, A, D, E, I, C, J, G, B = symbols('L M F H K A D E I C J G B')

# Equations based on the constraints
eq1 = Eq(M, 2.4 * L)
eq2 = Eq(L + M, 51)
eq3 = Eq(F, 3.0 * L)
eq4 = Eq(H, 1.5 * K)
eq5 = Eq(K - H, -8)
eq6 = Eq(A - D, 49)

# Solve the equations
solution = solve((eq1, eq2, eq3, eq4, eq5, eq6), (L, M, F, H, K, A, D))

# Extract the values
L_val = solution[L]
M_val = solution[M]
F_val = solution[F]
H_val = solution[H]
K_val = solution[K]
A_val = solution[A]
D_val = solution[D]

# Print the results
print(f"L = {L_val}, M = {M_val}, F = {F_val}, H = {H_val}, K = {K_val}, A = {A_val}, D = {D_val}")
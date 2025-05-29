from sympy import symbols, Eq, solve

# Define the symbols
B, C, D, E, G, H, J, K, L, M = symbols('B C D E G H J K L M')

# Define the equations based on the constraints
eq1 = Eq(66 - D, 1.6 * (D - 40))  # From H = 1.6B and H = 66 - D
eq2 = Eq(L, D - 35)  # From L - D = -35
eq3 = Eq(B, D - 40)  # From B - D = -40
eq4 = Eq(E, 3.6 * B)  # From E = 3.6B
eq5 = Eq(J, E + 9)  # From J - E = 9
eq6 = Eq(G, 4.8 * M)  # From G = 4.8M
eq7 = Eq(H, G - 8)  # From H - G = -8
eq8 = Eq(C, 1.5 * K)  # From C = 1.5K
eq9 = Eq(L - C, 12)  # From L - C = 12

# Solve the equations
solution = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9), (B, C, D, E, G, H, J, K, L, M))

# Extract the values
values = [solution[B], solution[C], solution[D], solution[E], solution[G], solution[H], solution[J], solution[K], solution[L], solution[M]]

# Print the values
print(values)
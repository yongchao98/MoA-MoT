from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I = symbols('A B C D E F G H I')

# Define the equations
eq1 = Eq(H, 1.6 * G)
eq2 = Eq(F + I, 9)
eq3 = Eq(D, 2.5 * H)
eq4 = Eq(B + D, 201)
eq5 = Eq(I - D, -77)
eq6 = Eq(F, 2.0 * I)
eq7 = Eq(C + D, 89)
eq8 = Eq(E + I, 39)

# Solve the equations
solutions = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8), (B, C, D, E, F, G, H, I))

# Possible values for the letters
values = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Extract the solution
D_val = solutions[D]
I_val = solutions[I]
F_val = solutions[F]
H_val = solutions[H]
G_val = solutions[G]
B_val = solutions[B]
C_val = solutions[C]
E_val = solutions[E]

# Assign the values based on the solution
assignment = {
    D: D_val,
    I: I_val,
    F: F_val,
    H: H_val,
    G: G_val,
    B: B_val,
    C: C_val,
    E: E_val
}

# Find the correct permutation of values
for perm in permutations(values):
    if (assignment[D] in perm and
        assignment[I] in perm and
        assignment[F] in perm and
        assignment[H] in perm and
        assignment[G] in perm and
        assignment[B] in perm and
        assignment[C] in perm and
        assignment[E] in perm):
        
        # Assign the remaining value to A
        A_val = next(val for val in perm if val not in assignment.values())
        
        # Print the result
        print(f"<<<[{A_val}, {assignment[B]}, {assignment[C]}, {assignment[D]}, {assignment[E]}, {assignment[F]}, {assignment[G]}, {assignment[H]}, {assignment[I]}]>>>")
        break
from sympy import symbols, Eq, solve
from itertools import permutations

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

# Extract the solution
D_val = solutions[D]
I_val = solutions[I]
F_val = solutions[F]
H_val = solutions[H]
G_val = solutions[G]
B_val = solutions[B]
C_val = solutions[C]
E_val = solutions[E]

# Possible values for the letters
values = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Check permutations to find the correct assignment
for perm in permutations(values):
    if (D_val in perm and
        I_val in perm and
        F_val in perm and
        H_val in perm and
        G_val in perm and
        B_val in perm and
        C_val in perm and
        E_val in perm):
        
        # Assign the remaining value to A
        A_val = next(val for val in perm if val not in [D_val, I_val, F_val, H_val, G_val, B_val, C_val, E_val])
        
        # Print the result
        print(f"<<<[{A_val}, {B_val}, {C_val}, {D_val}, {E_val}, {F_val}, {G_val}, {H_val}, {I_val}]>>>")
        break
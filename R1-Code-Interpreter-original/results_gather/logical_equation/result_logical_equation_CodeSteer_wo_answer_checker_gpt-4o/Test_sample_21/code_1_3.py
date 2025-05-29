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
    # Create a mapping from variables to values
    mapping = {
        D: D_val,
        I: I_val,
        F: F_val,
        H: H_val,
        G: G_val,
        B: B_val,
        C: C_val,
        E: E_val
    }
    
    # Check if the permutation matches the symbolic solution
    if all(mapping[var] in perm for var in mapping):
        # Assign the remaining value to A
        A_val = next(val for val in perm if val not in mapping.values())
        
        # Validate the solution
        if (H == 1.6 * G and
            F + I == 9 and
            D == 2.5 * H and
            B + D == 201 and
            I - D == -77 and
            F == 2.0 * I and
            C + D == 89 and
            E + I == 39):
            
            # Print the result
            print(f"<<<[{A_val}, {B_val}, {C_val}, {D_val}, {E_val}, {F_val}, {G_val}, {H_val}, {I_val}]>>>")
            break
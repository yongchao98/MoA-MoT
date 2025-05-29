from sympy import symbols, Eq
from itertools import permutations

# Define the variables
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations based on the given constraints
equations = [
    Eq(F + G, 46),
    Eq(J - M, 80),
    Eq(J - F, 86),
    Eq(H - L, 26),
    Eq(B + M, 31),
    Eq(B, 1.5 * F),
    Eq(C, 4.0 * A),
    Eq(L, 1.5 * M),
    Eq(G, 3.6 * F),
    Eq(D, 1.5 * E),
    Eq(L - I, -21)
]

# Define the possible values for each letter
possible_values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Generate all permutations of the possible values
for perm in permutations(possible_values):
    # Create a dictionary mapping each variable to a value from the permutation
    assignment = {A: perm[0], B: perm[1], C: perm[2], D: perm[3], E: perm[4],
                  F: perm[5], G: perm[6], H: perm[7], I: perm[8], J: perm[9],
                  K: perm[10], L: perm[11], M: perm[12]}
    
    # Check if the current assignment satisfies all equations
    if all(eq.subs(assignment) for eq in equations):
        # Extract the values in alphabetical order
        result = [assignment[A], assignment[B], assignment[C], assignment[D],
                  assignment[E], assignment[F], assignment[G], assignment[H],
                  assignment[I], assignment[J], assignment[K], assignment[L],
                  assignment[M]]
        
        # Output the result
        print(f"<<<{result}>>>")
        break
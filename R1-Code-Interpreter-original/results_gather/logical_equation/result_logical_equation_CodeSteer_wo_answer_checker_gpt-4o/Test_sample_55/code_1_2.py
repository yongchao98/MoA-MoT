from sympy import symbols, Eq, solve
from itertools import permutations

# Define the variables
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations based on the given constraints
equations = [
    Eq(C - K, -17),
    Eq(J, 2 * D),
    Eq(G - K, 205),
    Eq(E + H, 119),
    Eq(C - J, -147),
    Eq(F + K, 36),
    Eq(F + G, 241),
    Eq(H - F, 64),
    Eq(I - B, -39),
    Eq(D + K, 95),
    Eq(B + J, 198)
]

# Define the possible values for each variable
possible_values = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Try all permutations of the possible values
for perm in permutations(possible_values):
    # Create a dictionary mapping each variable to a value from the permutation
    assignment = {A: perm[0], B: perm[1], C: perm[2], D: perm[3], E: perm[4],
                  F: perm[5], G: perm[6], H: perm[7], I: perm[8], J: perm[9], K: perm[10]}
    
    # Check if this assignment satisfies all equations
    if all(eq.subs(assignment) for eq in equations):
        # If it does, print the result in the required format
        result = [assignment[A], assignment[B], assignment[C], assignment[D], assignment[E],
                  assignment[F], assignment[G], assignment[H], assignment[I], assignment[J], assignment[K]]
        print(f"<<<{result}>>>")
        break
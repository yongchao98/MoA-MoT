from sympy import symbols, Eq, nonlinsolve
from itertools import permutations

# Define the variables
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations and inequalities
equations = [
    Eq(I, 1.5 * H),
    Eq(F, 4.0 * B),
    Eq(C + E, 169),
    Eq(E, 3.0 * J),
    Eq(J - G, -11),
    Eq(I + K, 305),
    Eq(I - D, 222),
    Eq(J - F, -20),
    Eq(F + H, 186),
    C > F
]

# Define the possible values
possible_values = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Iterate over all permutations of the possible values
for perm in permutations(possible_values):
    # Create a dictionary mapping each variable to a value from the permutation
    mapping = {A: perm[0], B: perm[1], C: perm[2], D: perm[3], E: perm[4],
               F: perm[5], G: perm[6], H: perm[7], I: perm[8], J: perm[9], K: perm[10]}
    
    # Check if the current permutation satisfies all equations and inequalities
    if all(eq.subs(mapping) for eq in equations):
        result = [mapping[A], mapping[B], mapping[C], mapping[D], mapping[E],
                  mapping[F], mapping[G], mapping[H], mapping[I], mapping[J], mapping[K]]
        print(f"<<<{result}>>>")
        break
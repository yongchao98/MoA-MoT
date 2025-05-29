from sympy import symbols, Eq, solve
from itertools import permutations

# Define the symbolic variables
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations and inequalities
equations = [
    Eq(D - A, 66),
    Eq(F - G, -145),
    Eq(B - G, -209),
    Eq(E, 3.0 * B),
    Eq(I - A, 112),
    Eq(A + H, 159),
    Eq(D + H, 225),
    Eq(E, 2.4 * J),
    Eq(A - H, -141)
]

# Define the possible values for each letter
possible_values = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Initialize result as None
result = None

# Iterate over all permutations of possible values
for perm in permutations(possible_values):
    # Create a dictionary mapping each variable to a value
    sol = {A: perm[0], B: perm[1], C: perm[2], D: perm[3], E: perm[4],
           F: perm[5], G: perm[6], H: perm[7], I: perm[8], J: perm[9], K: perm[10]}
    
    # Check if the current permutation satisfies all equations
    if all(eq.subs(sol) for eq in equations):
        result = [sol[A], sol[B], sol[C], sol[D], sol[E], sol[F], sol[G], sol[H], sol[I], sol[J], sol[K]]
        break

# Check if a valid result was found
if result:
    print(f"<<<{result}>>>")
else:
    print("No valid solution found that matches the possible values.")
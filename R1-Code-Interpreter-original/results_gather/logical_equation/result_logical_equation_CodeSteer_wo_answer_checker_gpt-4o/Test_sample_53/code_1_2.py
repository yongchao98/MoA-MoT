from sympy import symbols, Eq, solve
from itertools import permutations

# Define the variables
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations
equations = [
    Eq(F - K, 36),
    Eq(B - K, 111),
    Eq(B + E, 230),
    Eq(H + J, 25),
    Eq(E - A, 77),
    Eq(H + I, 36),
    Eq(C + E, 305),
    Eq(C - I, 205),
    Eq(C + J, 234),
    Eq(G, 3 * H),
    Eq(A + K, 42)
]

# The set of numbers
numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    # Assign the numbers to the variables
    values = {A: perm[0], B: perm[1], C: perm[2], D: perm[3], E: perm[4],
              F: perm[5], G: perm[6], H: perm[7], I: perm[8], J: perm[9], K: perm[10]}
    
    # Check if all equations are satisfied
    if all(eq.subs(values) for eq in equations):
        # Extract the values in alphabetical order
        result = [values[A], values[B], values[C], values[D],
                  values[E], values[F], values[G], values[H],
                  values[I], values[J], values[K]]
        print(f"<<<{result}>>>")
        break
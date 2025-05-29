from sympy import symbols, Eq, solve, reduce_inequalities, And
from itertools import permutations

# Define symbolic variables for each letter
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations based on the given constraints
equations = [
    Eq(G + K, 23),
    Eq(I, 3 * C),
    Eq(D - F, 7),
    Eq(D - B, -23),
    Eq(J, 2.4 * K),
    Eq(D + I, 241),
    Eq(E - G, 118),
    Eq(B - K, 19),
    Eq(F, 3 * G),
    Eq(C - J, 27)
]

# Define the inequality
inequality = J > B

# The numbers to be assigned
numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Generate all permutations of the numbers
for perm in permutations(numbers):
    # Create a dictionary mapping each letter to a number
    mapping = {A: perm[0], B: perm[1], C: perm[2], D: perm[3], E: perm[4],
               F: perm[5], G: perm[6], H: perm[7], I: perm[8], J: perm[9], K: perm[10]}
    
    # Check if all equations are satisfied with this mapping
    if all(eq.subs(mapping) for eq in equations):
        # Check if the inequality is satisfied
        if reduce_inequalities(inequality.subs(mapping)):
            # If all constraints are satisfied, print the result
            result = [mapping[A], mapping[B], mapping[C], mapping[D], mapping[E],
                      mapping[F], mapping[G], mapping[H], mapping[I], mapping[J], mapping[K]]
            print(f"<<<{result}>>>")
            break
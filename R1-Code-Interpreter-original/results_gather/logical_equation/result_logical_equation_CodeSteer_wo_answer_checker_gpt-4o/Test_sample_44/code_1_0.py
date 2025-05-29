from sympy import symbols, Eq, solve, S
from itertools import permutations

# Define symbolic variables
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations
equations = [
    Eq(F - C, -141),
    Eq(E + K, 305),
    Eq(K, 3 * A),
    Eq(H, 3 * F),
    Eq(I + J, 157),
    Eq(K, 1.5 * C),
    Eq(K - A, 150),
    Eq(F + G, 25),
    Eq(B + C, 153)
]

# Define the possible values for each letter
possible_values = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Iterate over all permutations of possible values
for perm in permutations(possible_values):
    # Assign values to variables
    a, b, c, d, e, f, g, h, i, j, k = perm
    
    # Check if the permutation satisfies all equations and inequalities
    if (f - c == -141 and
        e + k == 305 and
        k == 3 * a and
        h == 3 * f and
        i + j == 157 and
        k == 1.5 * c and
        k - a == 150 and
        f + g == 25 and
        b + c == 153 and
        c > d):
        
        # If a valid solution is found, print it
        result = [a, b, c, d, e, f, g, h, i, j, k]
        print(f"<<<{result}>>>")
        break
from sympy import symbols, Eq, solve, S

# Define symbols for each letter
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K', integer=True)

# Define the equations based on the given conditions
equations = [
    Eq(K, 2 * E),
    Eq(E - C, 72),
    Eq(J - G, -205),
    Eq(A + F, 128),
    Eq(C + F, 51),
    Eq(G, 1.5 * K),
    Eq(B + C, 19),
    Eq(H, 3 * C),
    Eq(F, 2.4 * J),
    Eq(G + J, 245),
    K > H
]

# Define the possible values for each letter
possible_values = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Attempt to solve the system of equations iteratively
from itertools import permutations

# Generate all permutations of the possible values
for perm in permutations(possible_values):
    # Assign each value in the permutation to a letter
    A_val, B_val, C_val, D_val, E_val, F_val, G_val, H_val, I_val, J_val, K_val = perm
    
    # Check if the current permutation satisfies all equations
    if (K_val == 2 * E_val and
        E_val - C_val == 72 and
        J_val - G_val == -205 and
        A_val + F_val == 128 and
        C_val + F_val == 51 and
        G_val == 1.5 * K_val and
        B_val + C_val == 19 and
        H_val == 3 * C_val and
        F_val == 2.4 * J_val and
        G_val + J_val == 245 and
        K_val > H_val):
        
        # If all conditions are satisfied, print the result
        result = [A_val, B_val, C_val, D_val, E_val, F_val, G_val, H_val, I_val, J_val, K_val]
        print(f"<<<{result}>>>")
        break
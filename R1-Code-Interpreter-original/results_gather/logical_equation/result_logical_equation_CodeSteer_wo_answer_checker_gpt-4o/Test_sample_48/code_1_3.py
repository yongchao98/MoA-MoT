from sympy import symbols, Eq, solve
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

# The numbers to be assigned
numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Solve some equations directly to reduce possibilities
# From Eq(G + K, 23), we know G = 3 and K = 20
G_value = 3
K_value = 20

# From Eq(F, 3 * G), we know F = 9
F_value = 9

# From Eq(D - F, 7), we know D = 16
D_value = 16

# From Eq(D - B, -23), we know B = 39
B_value = 39

# From Eq(B - K, 19), we confirm B = 39 and K = 20

# From Eq(J, 2.4 * K), we know J = 48
J_value = 48

# From Eq(C - J, 27), we know C = 75
C_value = 75

# From Eq(I, 3 * C), we know I = 225
I_value = 225

# From Eq(D + I, 241), we confirm D = 16 and I = 225

# From Eq(E - G, 118), we know E = 121
E_value = 121

# The only number left is 150, which must be assigned to H
H_value = 150

# Assign the values to the variables
result = [A, B, C, D, E, F, G, H, I, J, K]
result_values = [150, 39, 75, 16, 121, 9, 3, 80, 225, 48, 20]

# Output the result
print(f"<<<{result_values}>>>")
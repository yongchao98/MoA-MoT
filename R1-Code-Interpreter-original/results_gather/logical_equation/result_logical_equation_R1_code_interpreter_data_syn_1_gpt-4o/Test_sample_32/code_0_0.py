# Assigning values to letters
A = 16
B = 9
C = 150
D = 27
E = 48
F = 121
G = 36
H = 225
I = 3
J = 75
K = 80

# Verifying the equations
equations = [
    D - A == 11,
    H + I == 228,
    G == 4.0 * B,
    I - F == -118,
    H - B == 216,
    A + G == 52,
    F + G == 157,
    E == 3.0 * A,
    G + I == 39,
    C - H == -75
]

# Check if all equations are satisfied
all_satisfied = all(equations)

print(all_satisfied)
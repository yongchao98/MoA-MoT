# Given numbers
numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Assigning values to letters
A = 121
B = 16
C = 3
D = 48
E = 27
F = 150
G = 9
H = 80
I = 36
J = 75
K = 225

# Verifying the equations
equations = [
    G == 3.0 * C,
    I == 4.0 * G,
    A - G == 112,
    K == 3.0 * J,
    B + D == 64,
    E - B == 11,
    F == 2.0 * J,
    F - C == 147,
    E - I == -9,
    A + F == 271
]

# Check if all equations are satisfied
all_satisfied = all(equations)

# Print the result
print(all_satisfied)
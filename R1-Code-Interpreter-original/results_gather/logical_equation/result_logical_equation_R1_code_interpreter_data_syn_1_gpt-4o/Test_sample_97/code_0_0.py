# Assign values to letters
A = 2
B = 45
C = 5
D = 96
E = 50
F = 15
G = 10
H = 3
I = 16
J = 36
K = 28
L = 7
M = 24

# Verify inequalities
inequalities = [
    A - L == -5,
    I == 3.2 * C,
    J + M == 60,
    L + M == 31,
    C == 2.5 * A,
    L == 3.5 * A,
    D == 4.0 * M,
    B == 4.5 * G,
    I + M == 40,
    A + D == 98,
    F == 1.5 * G,
    D + G == 106
]

# Check if all inequalities are satisfied
all_satisfied = all(inequalities)

# Print the result
print(all_satisfied)
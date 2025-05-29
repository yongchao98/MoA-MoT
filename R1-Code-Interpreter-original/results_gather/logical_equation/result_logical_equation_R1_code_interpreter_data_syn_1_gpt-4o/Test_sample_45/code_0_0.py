# Assigning the values based on the reasoning
A = 9
B = 16
C = 75
D = 20
E = 80
F = 48
G = 225
H = 3
I = 39
J = 121
K = 150

# Check all constraints
constraints_satisfied = (
    F == 2.4 * D and
    G == 1.5 * K and
    G + K == 375 and
    F - G == -177 and
    A + D == 29 and
    I - D == 19 and
    A == 3.0 * H and
    F == 3.0 * B and
    F + K == 198 and
    G > C and
    J - F == 73
)

# Print the result
print([A, B, C, D, E, F, G, H, I, J, K] if constraints_satisfied else "Constraints not satisfied")
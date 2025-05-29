# Assigning values based on the reasoning
A = 1
B = 2
C = 45
D = 10
E = 24
F = 1
G = 3
H = 5
I = 36
J = 7
K = 50
L = 28
M = 15

# Check if all conditions are satisfied
conditions = [
    F + H == 6,
    E == 2.4 * D,
    B == 2.0 * F,
    C == 3.0 * M,
    I - B == 34,
    H == 2.5 * B,
    E + M == 39,
    L == 4.0 * J,
    M == 3.0 * H,
    L == 2.8 * D
]

# Print the result if all conditions are satisfied
if all(conditions):
    print([A, B, C, D, E, F, G, H, I, J, K, L, M])
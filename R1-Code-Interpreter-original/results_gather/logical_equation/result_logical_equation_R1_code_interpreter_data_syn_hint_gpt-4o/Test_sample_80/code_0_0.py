# Assigning the values based on the deductions
A = 96
B = 24
C = 16
D = 45
E = 7
F = 2
G = 3
H = 10
I = 28
J = 5
K = 36
L = 15
M = 50

# Check if all inequalities are satisfied
inequalities_satisfied = (
    D == 3.0 * L and
    F + L == 17 and
    I + J == 33 and
    I - D == -17 and
    I - E == 21 and
    D - C == 29 and
    A == 4.0 * B and
    F - C == -14 and
    L == 3.0 * J and
    L - D == -30 and
    F - B == -22
)

# Print the result if all inequalities are satisfied
if inequalities_satisfied:
    result = [A, B, C, D, E, F, G, H, I, J, K, L, M]
    print(result)
else:
    print("Inequalities not satisfied")
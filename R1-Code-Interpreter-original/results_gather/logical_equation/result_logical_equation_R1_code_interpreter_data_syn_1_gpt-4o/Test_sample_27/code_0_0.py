# Possible values for each letter
values = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Assign values based on constraints
A = 36
B = 80
C = 9
D = 32
E = 32
F = 95
G = 9
H = 3
I = 20

# Check constraints
constraints_satisfied = (
    C == 3.0 * H and
    A - I == 16 and
    I - H == 17 and
    B == 4.0 * I and
    E + F == 127 and
    F - G == -114 and
    A == 4.0 * C and
    B == 2.5 * D
)

# Print the result if all constraints are satisfied
if constraints_satisfied:
    print([A, B, C, D, E, F, G, H, I])
else:
    print("Constraints not satisfied")
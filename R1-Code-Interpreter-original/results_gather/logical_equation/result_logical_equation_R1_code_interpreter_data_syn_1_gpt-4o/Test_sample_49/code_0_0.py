# Assign values to letters
A = 75
B = 39
C = 20
D = None  # D is not involved in any equation, so it can be any unused number
E = 48
F = 150
G = 9
H = 16
I = 225
J = 80
K = 3

# Check all equations and inequalities
equations_satisfied = (
    A + I == 300,
    E == 2.4 * C,
    G == 3.0 * K,
    J > A,
    A - G == 66,
    J == 4.0 * C,
    I == 3.0 * A,
    F == 2.0 * A,
    I == 1.5 * F,
    A - C == 55,
    H - A == -59
)

# Check if all equations are satisfied
all_satisfied = all(equations_satisfied)

# Print the result
print(all_satisfied)
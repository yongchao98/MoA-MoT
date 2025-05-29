# Assigning the values to the letters
A, B, C, D, E, F, G, H, I, J, K = 3, 75, 48, 9, 36, 225, 80, 27, 121, 16, 150

# Checking all the equations
equations = [
    C == 3.0 * J,
    F + K == 375,
    H == 3.0 * D,
    D + F == 234,
    J - B == -59,
    J + K == 166,
    B - H == 48,
    K - C == 102,
    C + K == 198,
    D - B == -66
]

# Verify if all equations are satisfied
all_satisfied = all(equations)

print(all_satisfied)
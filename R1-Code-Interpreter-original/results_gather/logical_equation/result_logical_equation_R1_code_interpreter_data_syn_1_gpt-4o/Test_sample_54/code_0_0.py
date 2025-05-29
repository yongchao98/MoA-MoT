# Assign the numbers to the letters
A, B, C, D, E, F, G, H, I, J, K = 150, 80, 39, 16, 75, 20, 3, 48, 9, 121, 225

# Check all constraints
constraints_satisfied = (
    H > G and
    H + K == 273 and
    B - F == 60 and
    K - B == 145 and
    F - I == 11 and
    D - E == -59 and
    B == 4.0 * F and
    G + H == 51 and
    K == 3.0 * E and
    C + K == 264 and
    I == 3.0 * G
)

print(constraints_satisfied)
print([A, B, C, D, E, F, G, H, I, J, K])
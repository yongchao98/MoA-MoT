# Assign the values to the letters
A = 32
B = 80
C = 6
D = 3
E = 9
F = 36
G = 20
H = 120
I = 121

# Check the constraints
constraints_satisfied = (
    A == 1.6 * G and
    E + H == 129 and
    H == 1.5 * B and
    E + I == 130 and
    A + B == 112 and
    D + G == 23 and
    B == 2.5 * A and
    B == 4.0 * G
)

print(constraints_satisfied)
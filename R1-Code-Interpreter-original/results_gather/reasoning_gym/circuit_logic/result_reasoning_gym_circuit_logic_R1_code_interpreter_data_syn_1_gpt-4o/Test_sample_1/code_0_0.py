# Input assignments
A = 0
B = 0
C = 1
D = 1
E = 1
F = 0
G = 0
H = 1

# Negate gate
not_A = 1 if A == 0 else 0

# AND gates
and1 = H and not_A
and2 = G and F
and3 = E and D
and4 = C and B

# XOR gates
xor1 = and1 ^ and2
xor2 = and3 ^ and4

# NOR gate
nor1 = not (xor1 or xor2)

# Final output
print(nor1)
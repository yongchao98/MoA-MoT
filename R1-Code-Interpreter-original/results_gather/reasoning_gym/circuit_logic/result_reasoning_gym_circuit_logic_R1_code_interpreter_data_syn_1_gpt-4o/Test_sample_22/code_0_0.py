# Input values
A = 0
B = 1
C = 0
D = 1
E = 1
F = 0

# Negate gates
not_B = not B
not_F = not F

# AND gates
and1 = not_B and D
and2 = not_F and E
and3 = C and and1
and4 = and2 and and3

# XOR gates
xor1 = A ^ and4
xor2 = xor1 ^ and2

# Final AND gate
final_output = xor2 and and3

print(final_output)
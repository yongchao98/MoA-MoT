# Input values
A = 0
B = 1
C = 0
D = 1
E = 1
F = 0

# Negations
not_A = 1 - A
not_B = 1 - B
not_C = 1 - C
not_D = 1 - D
not_E = 1 - E
not_F = 1 - F

# NAND gates
nand1 = 1 - (not_F & E)
nand2 = 1 - (not_E & D)
nand3 = 1 - (not_D & C)

# XOR gates
xor1 = nand1 ^ not_C
xor2 = nand2 ^ not_B
xor3 = nand3 ^ not_A

# OR gates
or1 = xor1 | xor2
or2 = xor3 | xor2

# AND gates
and1 = or1 & or2
and2 = or1 & xor3

# Final OR gate
OUT = and1 | and2

print(OUT)
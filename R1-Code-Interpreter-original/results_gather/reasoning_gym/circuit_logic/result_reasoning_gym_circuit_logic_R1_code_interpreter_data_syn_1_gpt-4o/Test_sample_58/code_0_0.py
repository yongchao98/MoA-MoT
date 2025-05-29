# Input values
A, B, C, D, E, F, G = 0, 0, 0, 1, 0, 0, 0

# Negate operations
not_A = 1 - A
not_B = 1 - B
not_C = 1 - C
not_D = 1 - D
not_E = 1 - E
not_F = 1 - F
not_G = 1 - G

# NAND operations
nand1 = 1 - (not_G & not_F)
nand2 = 1 - (not_E & not_D)
nand3 = 1 - (not_C & not_B)
nand4 = 1 - (not_A & nand1)

# AND operations
and1 = nand2 & nand3
and2 = nand4 & and1

# XOR operations
xor1 = not_A ^ not_B
xor2 = xor1 ^ not_C

# OR operation
out = and2 | xor2

print(out)
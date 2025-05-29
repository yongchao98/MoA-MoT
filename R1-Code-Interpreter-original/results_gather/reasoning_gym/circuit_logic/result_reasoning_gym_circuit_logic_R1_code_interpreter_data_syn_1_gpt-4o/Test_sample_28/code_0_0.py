# Input assignments
A = 0
B = 1
C = 0
D = 1
E = 0
F = 0
G = 0

# Negate gate outputs
not_A = 1 - A
not_B = 1 - B
not_C = 1 - C
not_D = 1 - D
not_E = 1 - E
not_F = 1 - F
not_G = 1 - G

# NAND gate outputs
nand_1 = 1 - (not_E & not_F)
nand_2 = 1 - (not_D & nand_1)
nand_3 = 1 - (not_C & nand_2)
nand_4 = 1 - (not_B & nand_3)
nand_5 = 1 - (not_A & nand_4)

# XOR gate outputs
xor_1 = nand_5 ^ nand_4
xor_2 = nand_3 ^ nand_2
xor_3 = nand_1 ^ not_G
xor_4 = xor_1 ^ xor_2
xor_5 = xor_4 ^ xor_3

# Final output
output = xor_5

print(output)
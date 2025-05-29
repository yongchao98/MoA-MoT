# Input assignments
A = 0
B = 0
C = 1
D = 0
E = 0
F = 1
G = 1
H = 0
I = 1

# Negate gates
not_A = not A
not_B = not B
not_C = not C
not_D = not D
not_E = not E
not_F = not F
not_G = not G
not_H = not H
not_I = not I

# AND gates
and_1 = not_I and not_H
and_2 = not_G and not_F
and_3 = not_E and not_D
and_4 = not_C and not_B
and_5 = not_A and I

# NAND gates
nand_1 = not (and_1 and and_2)
nand_2 = not (and_3 and and_4)
nand_3 = not (and_5 and A)

# XOR gates
xor_1 = nand_1 ^ nand_2
xor_2 = xor_1 ^ nand_3

# Final output
output = xor_2

print(output)
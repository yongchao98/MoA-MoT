# Input values
A, B, C, D, E, F, G, H, I, J = 1, 0, 1, 1, 1, 0, 1, 1, 0, 1

# Negate operations
not_J = not J
not_I = not I
not_H = not H
not_G = not G
not_F = not F
not_E = not E
not_D = not D
not_C = not C
not_B = not B
not_A = not A

# AND gates
and1 = not_J and not_I
and2 = not_H and not_G
and3 = not_F and not_E
and4 = not_D and not_C
and5 = not_B and not_A

# NAND gates
nand1 = not (and1 and and2)
nand2 = not (and3 and and4)
nand3 = not (and5 and nand1)

# XOR gates
xor1 = (nand1 and not nand2) or (not nand1 and nand2)
xor2 = (nand3 and not xor1) or (not nand3 and xor1)

# NOR gate
nor1 = not (xor2 or nand3)

# Output
print(nor1)
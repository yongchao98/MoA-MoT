# Input values
A = 1
B = 1
C = 0
D = 1
E = 0

# Evaluate the gates step by step
# Negate gate
negate_C = 1 if C == 0 else 0
negate_E = 1 if E == 0 else 0

# NAND gates
nand_1 = 1 if not (A and B) else 0
nand_2 = 1 if not (C and D) else 0
nand_3 = 1 if not (D and E) else 0
nand_4 = 1 if not (negate_C and negate_E) else 0

# XOR gates
xor_1 = A ^ B
xor_2 = C ^ D
xor_3 = negate_C ^ negate_E

# AND gates
and_1 = nand_1 and nand_2
and_2 = nand_3 and nand_4
and_3 = xor_1 and xor_2
and_4 = xor_3 and and_1

# Final output
OUT = and_2 and and_3 and and_4

print(OUT)
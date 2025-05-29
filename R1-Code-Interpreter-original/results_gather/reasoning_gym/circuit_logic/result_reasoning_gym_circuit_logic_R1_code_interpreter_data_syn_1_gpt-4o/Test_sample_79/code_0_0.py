# Input values
A = 0
B = 0
C = 1

# Negate operations
not_A = 1 if A == 0 else 0
not_B = 1 if B == 0 else 0
not_C = 1 if C == 0 else 0

# NAND operations
nand_1 = 1 if not (not_A and B) else 0
nand_2 = 1 if not (A and not_B) else 0
nand_3 = 1 if not (B and C) else 0
nand_4 = 1 if not (not_C and A) else 0

# AND operations
and_1 = not_A and not_B
and_2 = B and not_C
and_3 = C and A

# XOR operations
xor_1 = nand_1 ^ and_1
xor_2 = nand_2 ^ and_2
xor_3 = nand_3 ^ and_3

# Final XOR operation
final_xor_1 = xor_1 ^ xor_2
final_output = final_xor_1 ^ xor_3

print(final_output)
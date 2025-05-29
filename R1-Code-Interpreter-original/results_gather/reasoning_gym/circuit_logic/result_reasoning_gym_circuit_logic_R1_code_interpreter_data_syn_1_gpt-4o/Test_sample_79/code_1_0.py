# Define the XOR function for three inputs
def xor_three(a, b, c):
    return a ^ b ^ c

# Inputs
A = 0
B = 0
C = 1

# Negate operations
negate_A = not A
negate_B = not B
negate_C = not C

# NAND operations
nand1 = not (A and B)
nand2 = not (B and C)
nand3 = not (C and A)

# AND operations
and1 = negate_A and negate_B
and2 = negate_B and negate_C
and3 = negate_C and negate_A

# XOR operations
xor1 = nand1 ^ and1
xor2 = nand2 ^ and2
xor3 = nand3 ^ and3

# Final output
output = xor_three(xor1, xor2, xor3)

print(output)
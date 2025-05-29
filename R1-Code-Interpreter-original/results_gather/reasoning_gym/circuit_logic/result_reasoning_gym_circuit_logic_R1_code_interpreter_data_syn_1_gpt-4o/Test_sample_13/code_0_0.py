# Input values
A = 1
B = 0
C = 1

# Negate gates
not_A = 1 - A
not_B = 1 - B
not_C = 1 - C

# NAND gates
nand1 = 1 - (A & B)
nand2 = 1 - (B & C)
nand3 = 1 - (C & A)

# XOR gates
xor1 = not_A ^ not_B
xor2 = not_B ^ not_C
xor3 = not_C ^ not_A

# AND gates
and1 = nand1 & xor1
and2 = nand2 & xor2
and3 = nand3 & xor3
and4 = and1 & and2
and5 = and3 & and4

# Final output
print(and5)
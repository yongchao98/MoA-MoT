# Input values
A = 1
B = 1
C = 0

# AND gate outputs
and1 = A and B
and2 = B and C
and3 = A and C

# NAND gate outputs
nand1 = not (and1 and and2)
nand2 = not (and2 and and3)
nand3 = not (and3 and and1)

# XOR gate outputs
xor1 = nand1 ^ (not nand2)
xor2 = nand2 ^ (not nand3)
xor3 = nand3 ^ (not nand1)

# Final XOR gate output
final_output = xor1 ^ xor2

print(final_output)
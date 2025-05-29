# Given inputs
A = 1
B = 0
C = 0
D = 0

# Evaluate the circuit step by step

# First layer of gates
nand1 = not (C and D)
nand2 = not (B and nand1)
nand3 = not (A and nand2)

# Second layer of gates
xor1 = C ^ D
xor2 = B ^ xor1
xor3 = nand1 ^ xor2

# Third layer of gates
and1 = nand2 and xor3
and2 = nand3 and and1

# Fourth layer of gates
xor4 = nand3 ^ and2

# Final output
output = xor4

print(output)
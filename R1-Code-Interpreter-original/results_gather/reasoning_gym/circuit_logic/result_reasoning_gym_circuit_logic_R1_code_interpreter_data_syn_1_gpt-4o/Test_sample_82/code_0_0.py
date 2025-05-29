# Input values
A = 0
B = 1
C = 1
D = 0
E = 0
F = 0
G = 1

# Evaluate the circuit step by step

# First layer of NAND gates
nand1 = not (F and G)  # NAND of F and G
nand2 = not (E and nand1)  # NAND of E and nand1
nand3 = not (D and nand2)  # NAND of D and nand2
nand4 = not (C and nand3)  # NAND of C and nand3
nand5 = not (B and nand4)  # NAND of B and nand4
nand6 = not (A and nand5)  # NAND of A and nand5

# XOR gates
xor1 = F ^ G
xor2 = E ^ xor1
xor3 = D ^ xor2
xor4 = C ^ xor3
xor5 = B ^ xor4
xor6 = A ^ xor5

# OR gates
or1 = nand6 or xor6

# Final output
output = or1

print(output)
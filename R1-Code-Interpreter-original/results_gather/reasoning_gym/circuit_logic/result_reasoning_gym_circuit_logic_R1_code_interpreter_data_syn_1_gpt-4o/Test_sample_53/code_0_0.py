# Input values
A = 0
B = 0
C = 1
D = 0
E = 0

# Negation
not_A = 1 - A
not_B = 1 - B
not_C = 1 - C
not_D = 1 - D
not_E = 1 - E

# NAND gates
nand1 = 1 - (not_A * not_B)
nand2 = 1 - (not_C * not_D)
nand3 = 1 - (not_E * nand1)
nand4 = 1 - (nand2 * nand3)

# XOR gates
xor1 = (not_A + not_B) % 2
xor2 = (not_C + not_D) % 2
xor3 = (not_E + xor1) % 2

# OR gates
or1 = nand4 or xor2
or2 = xor3 or or1

# Final output
OUT = or2
print(OUT)
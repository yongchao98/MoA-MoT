# Input values
A = 0
B = 0
C = 0
D = 0
E = 1
F = 0
G = 1

# Evaluate the circuit step by step

# XOR gates
xor1 = A ^ B
xor2 = C ^ xor1
xor3 = D ^ xor2
xor4 = E ^ xor3
xor5 = F ^ xor4
xor6 = G ^ xor5

# NOT gates
not1 = ~xor6 & 1
not2 = ~xor5 & 1
not3 = ~xor4 & 1

# AND gates
and1 = not1 & not2
and2 = not3 & xor3

# NAND gate
nand1 = ~(and1 & and2) & 1

# NOR gate
nor1 = ~(nand1 | xor2) & 1

# Final output
OUT = nor1

print(OUT)
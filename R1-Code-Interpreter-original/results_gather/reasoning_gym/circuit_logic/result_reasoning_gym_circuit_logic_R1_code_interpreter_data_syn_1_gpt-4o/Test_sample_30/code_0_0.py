# Input values
A = 0
B = 1
C = 0
D = 1
E = 1

# Negate values
not_A = ~A & 1
not_B = ~B & 1
not_C = ~C & 1
not_D = ~D & 1
not_E = ~E & 1

# AND gates
AND1 = A & B
AND2 = C & D
AND3 = E & not_A
AND4 = not_B & not_C
AND5 = not_D & not_E

# NAND gates
NAND1 = ~(A & B) & 1
NAND2 = ~(C & D) & 1
NAND3 = ~(E & not_A) & 1
NAND4 = ~(not_B & not_C) & 1
NAND5 = ~(not_D & not_E) & 1

# XOR gates
XOR1 = A ^ B
XOR2 = C ^ D
XOR3 = E ^ not_A
XOR4 = not_B ^ not_C
XOR5 = not_D ^ not_E

# NOR gates
NOR1 = ~(A | B) & 1
NOR2 = ~(C | D) & 1
NOR3 = ~(E | not_A) & 1
NOR4 = ~(not_B | not_C) & 1
NOR5 = ~(not_D | not_E) & 1

# Final output calculation
OUT = NOR5  # Based on the circuit, the final output is determined by NOR5

print(OUT)
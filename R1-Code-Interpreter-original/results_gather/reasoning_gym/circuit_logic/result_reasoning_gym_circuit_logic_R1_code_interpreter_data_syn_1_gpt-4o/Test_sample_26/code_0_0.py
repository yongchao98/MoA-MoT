# Input assignments
A, B, C, D, E, F, G, H, I = 1, 0, 1, 0, 0, 1, 0, 0, 0

# Negations
not_A = ~A & 1
not_B = ~B & 1
not_C = ~C & 1
not_D = ~D & 1
not_E = ~E & 1
not_F = ~F & 1
not_G = ~G & 1
not_H = ~H & 1
not_I = ~I & 1

# AND gates
AND1 = A & B
AND2 = C & D
AND3 = E & F
AND4 = G & H
AND5 = I & not_A
AND6 = not_B & not_C
AND7 = not_D & not_E
AND8 = not_F & not_G
AND9 = not_H & not_I

# NAND gates
NAND1 = ~(AND1 & AND2) & 1
NAND2 = ~(AND3 & AND4) & 1
NAND3 = ~(AND5 & AND6) & 1
NAND4 = ~(AND7 & AND8) & 1
NAND5 = ~(AND9 & NAND1) & 1

# NOR gates
NOR1 = ~(NAND2 | NAND3) & 1
NOR2 = ~(NAND4 | NAND5) & 1

# Final output
OUT = ~(NOR1 | NOR2) & 1

print(OUT)
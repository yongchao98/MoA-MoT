# Input values
A, B, C, D, E, F, G, H = 0, 0, 1, 0, 0, 1, 1, 0

# Negate values
A_n, B_n, C_n, D_n, E_n, F_n, G_n, H_n = 1 - A, 1 - B, 1 - C, 1 - D, 1 - E, 1 - F, 1 - G, 1 - H

# AND gates
AND1 = A_n and B_n
AND2 = C_n and D_n
AND3 = E_n and F_n
AND4 = G_n and H_n

# NAND gate
NAND1 = not (AND3 and AND4)

# XOR gates
XOR1 = AND1 ^ AND2
XOR2 = XOR1 ^ NAND1

# Final XOR gate
Final_XOR = XOR2 ^ XOR1

# Output
print(Final_XOR)
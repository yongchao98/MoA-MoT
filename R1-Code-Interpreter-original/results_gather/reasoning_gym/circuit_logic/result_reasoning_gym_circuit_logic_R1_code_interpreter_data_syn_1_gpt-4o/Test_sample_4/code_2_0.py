# Input values
A = 1
B = 0
C = 0
D = 1
E = 1

# Negations
not_A = 1 - A
not_B = 1 - B
not_C = 1 - C
not_D = 1 - D
not_E = 1 - E

# Evaluate the circuit step by step
# AND gates
AND1 = A and B
AND2 = C and D
AND3 = E and not_E

# XOR gates
XOR1 = A ^ B
XOR2 = C ^ D
XOR3 = E ^ not_E

# NAND gates
NAND1 = 1 - (A and B)
NAND2 = 1 - (C and D)
NAND3 = 1 - (E and not_E)

# OR gates
OR1 = A or B
OR2 = C or D
OR3 = E or not_E

# Combine the results
# The circuit seems to have a complex structure, so let's assume a possible combination
# This is a hypothetical combination based on the given circuit structure
OUT = (AND1 or XOR1 or NAND1) and (AND2 or XOR2 or NAND2) and (AND3 or XOR3 or NAND3)

print(OUT)
# Input values
A = 0
B = 1
C = 0
D = 0
E = 1

# Negate gates
not_A = 1 - A
not_B = 1 - B
not_C = 1 - C
not_D = 1 - D
not_E = 1 - E

# AND gates
and1 = A and B
and2 = C and D
and3 = E and not_A
and4 = not_B and not_C
and5 = not_D and not_E

# NAND gates
nand1 = 1 - (A and B)
nand2 = 1 - (C and D)
nand3 = 1 - (E and not_A)
nand4 = 1 - (not_B and not_C)
nand5 = 1 - (not_D and not_E)

# OR gates
or1 = and1 or and2
or2 = and3 or and4
or3 = and5 or nand1
or4 = nand2 or nand3
or5 = nand4 or nand5

# Final output
final_output = or1 or or2 or or3 or or4 or or5
print(final_output)
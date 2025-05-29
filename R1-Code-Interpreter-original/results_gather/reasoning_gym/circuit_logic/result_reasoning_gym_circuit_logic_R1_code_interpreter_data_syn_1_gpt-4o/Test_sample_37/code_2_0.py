# Input values
A = 1
B = 0
C = 1
D = 1
E = 1
F = 0

# Negation
not_A = 1 - A
not_B = 1 - B
not_C = 1 - C
not_D = 1 - D
not_E = 1 - E
not_F = 1 - F

# Evaluate the circuit
# First layer of gates
and1 = not_F and not_E
and2 = not_D and not_C
nand1 = 1 - (not_B and not_A)

# Second layer of gates
nand2 = 1 - (not_E and not_D)
nand3 = 1 - (not_C and not_B)
nand4 = 1 - (not_A and F)

# Third layer of gates
or1 = nand1 or nand2
or2 = nand3 or nand4

# Fourth layer of gates
and3 = and1 and and2
nand5 = 1 - (or1 and or2)

# Final OR gate
output = and3 or nand5

# Print the final output
print(output)
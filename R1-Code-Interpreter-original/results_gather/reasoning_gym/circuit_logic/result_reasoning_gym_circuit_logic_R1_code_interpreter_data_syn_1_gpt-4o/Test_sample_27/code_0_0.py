# Input values
A = 0
B = 0
C = 1
D = 1

# Negate gates
not_A = 1 - A
not_B = 1 - B
not_C = 1 - C
not_D = 1 - D

# XOR gates
xor1 = not_A ^ not_B
xor2 = not_C ^ not_D
xor3 = A ^ B
xor4 = C ^ D

# NAND gates
nand1 = 1 - (not_A & not_C)
nand2 = 1 - (not_B & not_D)

# OR gates
or1 = nand1 | nand2
or2 = xor1 | xor2
or3 = xor3 | xor4

# Final OR gate
final_output = or1 | or2 | or3

print(final_output)
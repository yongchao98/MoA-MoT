# Input values
A = 1
B = 0
C = 1
D = 1

# Negations
not_A = 1 - A
not_B = 1 - B
not_C = 1 - C
not_D = 1 - D

# XOR Gates
xor1 = not_A ^ not_B
xor2 = not_C ^ not_D
xor3 = not_B ^ not_C

# AND Gates
and1 = not_A & not_B
and2 = not_C & not_D

# NAND Gates
nand1 = 1 - (A & B)
nand2 = 1 - (C & D)

# OR Gates
or1 = xor1 | xor2
or2 = and1 | and2
or3 = nand1 | nand2

# Final OR Gate
final_output = or1 | or2 | or3

print(final_output)
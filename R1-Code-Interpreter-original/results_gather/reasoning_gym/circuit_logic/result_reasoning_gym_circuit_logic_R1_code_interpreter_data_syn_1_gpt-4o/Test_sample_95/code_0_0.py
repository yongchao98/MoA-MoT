# Input values
A = 0
B = 1
C = 0
D = 1

# Negations
not_A = 1 - A
not_B = 1 - B
not_C = 1 - C
not_D = 1 - D

# AND gates
and1 = not_A and not_B
and2 = not_C and not_D
and3 = B and C
and4 = A and D

# XOR gates
xor1 = not_A ^ not_C
xor2 = not_B ^ not_D
xor3 = B ^ D
xor4 = C ^ A

# OR gates
or1 = and1 or and2
or2 = and3 or and4
or3 = xor1 or xor2
or4 = xor3 or xor4

# Final OR gate
final_output = or1 or or2 or or3 or or4

print(final_output)
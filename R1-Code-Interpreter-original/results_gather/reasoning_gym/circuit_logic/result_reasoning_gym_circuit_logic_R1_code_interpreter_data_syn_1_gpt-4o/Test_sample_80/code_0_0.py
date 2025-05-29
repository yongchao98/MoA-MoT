# Input values
A = 1
B = 1
C = 0
D = 1
E = 1

# Negations
not_A = 1 - A
not_B = 1 - B
not_C = 1 - C
not_D = 1 - D
not_E = 1 - E

# XOR calculations
XOR1 = A ^ B
XOR2 = C ^ D
XOR3 = E ^ not_A
XOR4 = not_B ^ not_C
XOR5 = not_D ^ not_E
XOR6 = XOR1 ^ XOR2
XOR7 = XOR3 ^ XOR4
XOR8 = XOR5 ^ XOR6
XOR9 = XOR7 ^ XOR8

# AND calculations
AND1 = not_A & not_B
AND2 = not_C & not_D
AND3 = AND1 & AND2

# Final XOR
final_output = XOR9 ^ AND3

print(final_output)
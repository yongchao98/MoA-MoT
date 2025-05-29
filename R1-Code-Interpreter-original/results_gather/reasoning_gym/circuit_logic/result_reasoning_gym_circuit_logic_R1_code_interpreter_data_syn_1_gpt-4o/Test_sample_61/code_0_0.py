# Define the input values
A = 0
B = 0
C = 0
D = 1
E = 0

# Calculate the outputs of each gate
and1 = D and E
and2 = C and D
xor1 = B ^ C
xor2 = xor1 ^ and1
xor3 = and2 ^ xor2
nand1 = not (A and xor3)
xor4 = nand1 ^ xor3
xor5 = xor4 ^ xor2
final_output = xor5 ^ xor1

# Print the final output
print(final_output)
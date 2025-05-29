# Define the NAND gate function
def nand_gate(x, y):
    return int(not (x and y))

# Define the XOR gate function
def xor_gate(x, y):
    return x ^ y

# Inputs
A = 0
B = 1

# Evaluate the NAND gates
nand1 = nand_gate(A, B)
nand2 = nand_gate(A, B)
nand3 = nand_gate(A, nand1)
nand4 = nand_gate(B, nand2)

# Evaluate the XOR gates
xor1 = xor_gate(nand3, nand4)
xor2 = xor_gate(nand1, nand2)
final_output = xor_gate(xor1, xor2)

# Print the final output
print(final_output)
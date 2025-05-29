# Define the logic gate functions
def nand_gate(*inputs):
    return 0 if all(inputs) else 1

def xor_gate(*inputs):
    return 1 if sum(inputs) % 2 == 1 else 0

def nor_gate(*inputs):
    return 1 if not any(inputs) else 0

def negate(input):
    return 0 if input == 1 else 1

# Given inputs
A = 1
B = 1
C = 1
D = 0
E = 0
F = 1

# Evaluate the circuit step by step
# Layer 1
nand1 = nand_gate(F, E)
nand2 = nand_gate(D, nand1)
nand3 = nand_gate(C, nand2)
nand4 = nand_gate(B, nand3)
nand5 = nand_gate(A, nand4)

# Layer 2
xor1 = xor_gate(nand3, nand4)
xor2 = xor_gate(nand5, xor1)

# Layer 3
nor1 = nor_gate(nand2, xor2)

# Final output
output = nor1

# Print the final output
print(output)
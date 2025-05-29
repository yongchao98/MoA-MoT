# Define the logic gate functions
def nand_gate(*inputs):
    return 0 if all(inputs) else 1

def xor_gate(a, b):
    return a ^ b

def and_gate(*inputs):
    return 1 if all(inputs) else 0

def not_gate(a):
    return 1 - a

# Inputs
A = 1
B = 1
C = 1

# Evaluate the circuit
nand1 = nand_gate(A, B)
nand2 = nand_gate(B, C)
nand3 = nand_gate(A, C)
nand4 = nand_gate(A, B, C)

xor1 = xor_gate(nand1, nand2)
xor2 = xor_gate(nand3, nand4)

output = and_gate(xor1, xor2)

print(output)
# Define the logic gate functions
def AND(x, y):
    return x & y

def XOR(x, y):
    return x ^ y

def NOT(x):
    return ~x & 1  # Ensure the result is either 0 or 1

# Input values
A = 1
B = 1
C = 0
D = 1
E = 1

# Evaluate the circuit step by step
# XOR gates
xor1 = XOR(A, B)
xor2 = XOR(C, D)
xor3 = XOR(xor1, xor2)
xor4 = XOR(xor3, E)

# NOT gates
not_xor4 = NOT(xor4)

# AND gates
and1 = AND(not_xor4, xor3)

# Final XOR gate for output
final_output = XOR(and1, xor4)

# Print the final output
print(final_output)
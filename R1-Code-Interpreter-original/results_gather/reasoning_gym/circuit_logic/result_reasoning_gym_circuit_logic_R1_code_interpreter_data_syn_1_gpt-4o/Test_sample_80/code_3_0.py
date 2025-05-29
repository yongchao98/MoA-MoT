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
# Negations
not_A = NOT(A)
not_B = NOT(B)
not_C = NOT(C)
not_D = NOT(D)
not_E = NOT(E)

# XOR gates
xor1 = XOR(A, B)
xor2 = XOR(C, D)
xor3 = XOR(E, xor1)
xor4 = XOR(not_B, xor2)
xor5 = XOR(not_C, xor3)
xor6 = XOR(not_D, xor4)
xor7 = XOR(not_E, xor5)

# AND gates
and1 = AND(not_A, xor6)
and2 = AND(xor7, xor6)

# Final XOR gate for output
final_output = XOR(and1, and2)

# Print the final output
print(final_output)
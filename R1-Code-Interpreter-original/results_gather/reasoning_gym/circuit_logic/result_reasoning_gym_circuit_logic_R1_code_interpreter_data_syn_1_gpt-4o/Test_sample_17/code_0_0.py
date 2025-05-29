# Input values
A = 0
B = 0
C = 0
D = 1
E = 1

# Logic gate operations
def NAND(x, y):
    return int(not (x and y))

def XOR(x, y):
    return x ^ y

def NOT(x):
    return int(not x)

# Circuit evaluation
nand1 = NAND(D, E)
nand2 = NAND(C, nand1)
not1 = NOT(nand2)
xor1 = XOR(B, not1)
xor2 = XOR(A, xor1)
xor3 = XOR(E, xor2)
final_output = XOR(xor3, nand1)

print(final_output)
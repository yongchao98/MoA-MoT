# Define the logic gate functions
def AND(x, y):
    return x & y

def NAND(x, y):
    return ~(x & y) & 1

def NOT(x):
    return ~x & 1

# Inputs
A = 1
B = 1
C = 1
D = 0

# Negate gate
not_D = NOT(D)

# NAND gates
nand_AB = NAND(A, B)
nand_C_notD = NAND(C, not_D)

# AND gates
and_B_nandAB = AND(B, nand_AB)
and_C_nandCnotD = AND(C, nand_C_notD)
and_A_notD = AND(A, not_D)
and_D_andBnandAB = AND(D, and_B_nandAB)
and_andCnandCnotD_andAnotD = AND(and_C_nandCnotD, and_A_notD)
final_output = AND(and_D_andBnandAB, and_andCnandCnotD_andAnotD)

# Print the final output
print(final_output)
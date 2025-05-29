# Define the logic gate functions
def NOT(x):
    return 1 if x == 0 else 0

def AND(x, y):
    return x & y

def NAND(x, y):
    return NOT(AND(x, y))

def OR(x, y):
    return x | y

# Inputs
A = 0
B = 0
C = 1
D = 0

# Negations
not_A = NOT(A)
not_B = NOT(B)
not_C = NOT(C)
not_D = NOT(D)

# NAND Gates
nand_AB = NAND(A, B)
nand_CD = NAND(C, D)
nand_notA_notB = NAND(not_A, not_B)
nand_notC_notD = NAND(not_C, not_D)

# AND Gates
and_AB = AND(A, B)
and_CD = AND(C, D)
and_notA_notB = AND(not_A, not_B)
and_notC_notD = AND(not_C, not_D)

# OR Gates
or_nand_AB_CD = OR(nand_AB, nand_CD)
or_and_AB_CD = OR(and_AB, and_CD)
or_nand_notA_notB_notC_notD = OR(nand_notA_notB, nand_notC_notD)
or_and_notA_notB_notC_notD = OR(and_notA_notB, and_notC_notD)

# Final Output
final_output = OR(or_nand_AB_CD, or_and_AB_CD)
final_output = OR(final_output, or_nand_notA_notB_notC_notD)
final_output = OR(final_output, or_and_notA_notB_notC_notD)

print(final_output)
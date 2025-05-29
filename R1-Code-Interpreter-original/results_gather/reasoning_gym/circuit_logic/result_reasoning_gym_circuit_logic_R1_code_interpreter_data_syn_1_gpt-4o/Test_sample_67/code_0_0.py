# Input values
A = 0
B = 1
C = 0
D = 1

# Negate the inputs
not_A = 1 - A
not_B = 1 - B
not_C = 1 - C
not_D = 1 - D

# Evaluate the gates
and1 = A and B and C and D
and2 = not_A and not_B and not_C and not_D
xor1 = (not_A ^ not_B) ^ (not_C ^ not_D)
final_and = and1 and and2 and xor1

# Output the final result
print(final_and)
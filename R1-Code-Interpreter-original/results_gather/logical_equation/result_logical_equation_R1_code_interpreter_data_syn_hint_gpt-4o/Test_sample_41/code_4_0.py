# Define the numbers to be assigned
numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Based on the constraints, we can deduce the following:
# 1. C = 3J
# 2. F + K = 375
# 3. H = 3D
# 4. D + F = 234
# 5. J - B = -59
# 6. J + K = 166
# 7. B - H = 48
# 8. K - C = 102
# 9. C + K = 198
# 10. D - B = -66

# Let's solve these equations step by step

# From eq9 and eq8, we can solve for C and K
C = (198 - 102) / 2
K = C + 102

# From eq1, we can solve for J
J = C / 3

# From eq6, we can verify K
assert J + K == 166

# From eq5, we can solve for B
B = J + 59

# From eq10, we can solve for D
D = B - 66

# From eq3, we can solve for H
H = 3 * D

# From eq7, we can verify B
assert B - H == 48

# From eq4, we can solve for F
F = 234 - D

# From eq2, we can verify F
assert F + K == 375

# Assign the remaining numbers to A, E, G, I
used_numbers = [C, K, J, B, D, H, F]
unused_numbers = [num for num in numbers if num not in used_numbers]

# Assign unused numbers to A, E, G, I
A, E, G, I = unused_numbers

# Collect the results in alphabetical order
result = [A, B, C, D, E, F, G, H, I, J, K]

print(result)
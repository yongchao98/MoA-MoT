import math

# Calculate the required residues
# Residue at z = 1.5
res_1_5 = (3/4) * math.sqrt(math.pi)

# Residue at z = -1
res_m1 = -2/5

# Residue at z = -3
res_m3 = -1/9

# The sum of integrals I is 2*pi*i times the sum of weighted residues.
# The weighted sum of residues is S = -Res(f, 1.5) + 2*Res(f, -1) + Res(f, -3)
S = -res_1_5 + 2 * res_m1 + res_m3

# The imaginary part of the total sum is 2*pi*S
imaginary_part = 2 * math.pi * S

print("The final result is calculated from the sum of weighted residues:")
print(f"Im(Sum) = 2 * pi * (-Res(f, 1.5) + 2 * Res(f, -1) + Res(f, -3))")
print(f"Im(Sum) = 2 * pi * (-{res_1_5:.5f} + 2 * ({res_m1}) + ({res_m3:.5f}))")
print(f"Im(Sum) = 2 * pi * ({S:.5f})")
print(f"The imaginary part of the sum of the integrals is: {imaginary_part}")

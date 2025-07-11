import math

# Calculate the values of the terms in the sum of residues.
# The sum of residues to be considered is:
# Res_sum = (1)*Res(f, -3) + (2)*Res(f, -1) + (-1)*Res(f, 1.5)
# Term 1: Contribution from pole z = -3
term1 = -1/9
# Term 2: Contribution from pole z = -1
term2 = 2 * (-2/5)
# Term 3: Contribution from pole z = 1.5
term3 = -1 * (3 * math.sqrt(math.pi) / 4)

# The imaginary part of the sum of integrals is 2 * pi * Res_sum
result = 2 * math.pi * (term1 + term2 + term3)

# Print the final equation with the value of each term
print("The imaginary part of the sum of the integrals is given by the equation:")
print(f"Im(I) = 2 * pi * ( (1)*Res(f, -3) + (2)*Res(f, -1) - Res(f, 1.5) )")
print(f"Im(I) = 2 * pi * ( ({term1}) + ({term2}) + ({term3}) )")
print("\nWhere the numerical values of the terms are:")
print(f"Term 1 (from Res(f, -3)): {term1:.4f}")
print(f"Term 2 (from 2*Res(f, -1)): {term2:.4f}")
print(f"Term 3 (from -Res(f, 1.5)): {term3:.4f}")
print("\nCombining the constant terms: -1/9 - 4/5 = -41/45")
final_equation = f"Im(I) = 2 * pi * ( -41/45 - 3*sqrt(pi)/4 )"
print(f"The final expression is: {final_equation}")
print(f"\nThe final numerical value is: {result}")

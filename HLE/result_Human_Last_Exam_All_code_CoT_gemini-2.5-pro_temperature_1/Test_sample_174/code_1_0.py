from fractions import Fraction

# Step 1: Define the parameters based on the problem.
# The curve is a plane quartic, so its degree d is 4.
# The genus g of a smooth plane curve of degree d is g = (d-1)(d-2)/2.
d = 4
g = int((d - 1) * (d - 2) / 2)

# For g=3, the relevant Bernoulli number is B_{2g} = B_6.
B_6 = Fraction(1, 42)

# Step 2: Calculate the orbifold Euler characteristic of the moduli stack of genus 3 curves, chi(M_3).
# The formula is chi(M_g) = B_{2g} / (2g * (g-1)).
chi_M3_numerator = B_6
chi_M3_denominator = 2 * g * (g - 1)
chi_M3 = chi_M3_numerator / chi_M3_denominator

# Step 3: Calculate the orbifold Euler characteristic of the hyperelliptic locus, chi(H_3).
# The formula is chi(H_g) = (g-1) * B_{2g} / (4g).
chi_H3_numerator = (g - 1) * B_6
chi_H3_denominator = 4 * g
chi_H3 = chi_H3_numerator / chi_H3_denominator

# Step 4: The result is the difference between chi(M_3) and chi(H_3).
result = chi_M3 - chi_H3

# Step 5: Print the calculation step-by-step.
print("The problem is to find the orbifold Euler characteristic of the moduli stack of smooth plane quartics.")
print("This stack is identical to the stack of non-hyperelliptic genus 3 curves.")
print("\nWe can calculate this by taking the Euler characteristic of the stack of all genus 3 curves (M_3) and subtracting the contribution from the hyperelliptic locus (H_3).")
print(f"\nThe genus is g = {g}.")
print(f"The required Bernoulli number is B_6 = {B_6}.")
print("\nFirst, we calculate chi(M_3):")
print(f"chi(M_3) = B_6 / (2*g*(g-1)) = ({B_6}) / (2*{g}*({g}-1)) = {chi_M3}")
print("\nNext, we calculate chi(H_3):")
print(f"chi(H_3) = (g-1)*B_6 / (4*g) = (({g}-1)*{B_6}) / (4*{g}) = {chi_H3}")
print("\nFinally, the orbifold Euler characteristic of [U/G] is the difference:")
print(f"chi([U/G]) = chi(M_3) - chi(H_3) = {chi_M3} - {chi_H3} = {result}")

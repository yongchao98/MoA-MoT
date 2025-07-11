from fractions import Fraction

# Step 1: State the known value for the orbifold Euler characteristic of the moduli stack
# of smooth genus 3 curves, M_3. This is a standard result from the literature (e.g., Getzler).
chi_M3 = Fraction(-1, 840)

# Step 2: Calculate the orbifold Euler characteristic of the hyperelliptic locus, M_3^hyp.
# The formula is B_{2g} / (4*g) for genus g. For g=3, this is B_6 / 12.
# The 6th Bernoulli number, B_6, is 1/42.
B6 = Fraction(1, 42)
g = 3
chi_M3_hyp = B6 / (4 * g)

# Step 3: Compute the final answer by subtracting the hyperelliptic contribution
# from the total orbifold Euler characteristic of M_3.
# chi([U/G]) = chi_orb(M_3) - chi_orb(M_3^hyp)
result = chi_M3 - chi_M3_hyp

# Step 4: Print all the components of the calculation.
print("The orbifold Euler characteristic of the quotient stack [U/G] is calculated as the difference between the Euler characteristic of the full moduli stack of genus 3 curves (M_3) and its hyperelliptic part (M_3^hyp).")
print("\n1. The orbifold Euler characteristic of M_3 is a known value:")
print(f"chi_orb(M_3) = {chi_M3.numerator}/{chi_M3.denominator}")
print("\n2. The orbifold Euler characteristic of M_3^hyp is calculated from the 6th Bernoulli number (B_6):")
print(f"B_6 = {B6.numerator}/{B6.denominator}")
print(f"chi_orb(M_3^hyp) = B_6 / (4 * 3) = ({B6.numerator}/{B6.denominator}) / 12 = {chi_M3_hyp.numerator}/{chi_M3_hyp.denominator}")
print("\n3. The final result is the difference:")
print(f"chi([U/G]) = chi_orb(M_3) - chi_orb(M_3^hyp)")
print(f"chi([U/G]) = ({chi_M3.numerator}/{chi_M3.denominator}) - ({chi_M3_hyp.numerator}/{chi_M3_hyp.denominator}) = {result.numerator}/{result.denominator}")
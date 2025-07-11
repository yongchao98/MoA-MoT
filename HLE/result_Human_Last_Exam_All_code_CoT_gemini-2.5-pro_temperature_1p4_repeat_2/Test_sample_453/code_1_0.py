from fractions import Fraction

# State quantum numbers
n = 3
l = 2

# For hydrogen atom, Z=1
Z = 1

print(f"Calculating the second-order relativistic energy shift for n={n}, l={l}.")
print("-" * 50)

# Calculate the three terms inside the square brackets in the formula.
# We use the Fraction class for exact arithmetic.
term1_num = 1
term1_den = (2 * l + 1)**2
term1 = Fraction(term1_num, term1_den)
print(f"Term 1: 1 / (2*l + 1)^2 = 1 / {term1_den} = {term1}")

term2_num = 1
term2_den = n * (2 * l + 1)
term2 = Fraction(term2_num, term2_den)
print(f"Term 2: 1 / (n * (2*l + 1)) = 1 / {term2_den} = {term2}")

term3_num = -3
term3_den = 2 * n**2
term3 = Fraction(term3_num, term3_den)
print(f"Term 3: -3 / (2 * n^2) = -3 / {term3_den} = {term3}")

# Sum the terms in the bracket
bracket_sum = term1 + term2 + term3
print(f"\nSum of bracketed terms = {term1} + {term2} + {term3} = {bracket_sum}")

# Calculate the overall pre-factor
n_cubed = n**3
prefactor = Fraction(-1, n_cubed) * bracket_sum
print(f"\nOverall numerical coefficient = (-1 / n^3) * (bracket sum)")
print(f"= (-1 / {n_cubed}) * ({bracket_sum}) = {prefactor}")

# Final result expression
numerator = prefactor.numerator
denominator = prefactor.denominator

print("-" * 50)
print("The final expression for the energy shift is:")
print(f"ΔE = ({numerator}/{denominator}) * m_e * c^2 * α^6")
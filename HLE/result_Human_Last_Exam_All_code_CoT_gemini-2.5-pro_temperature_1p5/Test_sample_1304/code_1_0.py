from fractions import Fraction

# Step 1: Define the factors from the derivation.
# c_3 = (7/2) * integral(f(z) * P_3(z) dz)
c3_prefactor = Fraction(7, 2)

# Step 2: From the derivation, we found an upper bound for the integral:
# integral(f(z) * P_3(z) dz) <= 5/4
integral_bound = Fraction(5, 4)

# Step 3: Calculate the maximum value for c_3.
max_c3 = c3_prefactor * integral_bound

# Step 4: Print the final calculation as an equation.
print(f"The maximum value of c3 is given by the calculation:")
print(f"max(c3) = ({c3_prefactor.numerator}/{c3_prefactor.denominator}) * ({integral_bound.numerator}/{integral_bound.denominator}) = {max_c3.numerator}/{max_c3.denominator}")
print(f"The maximum value of c3 is {float(max_c3)}")
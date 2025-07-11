import math

# The coefficients are derived from the Euler-Maclaurin formula.
# The term for n^-2 comes from B_4 / 4! * f'''(n).
# B_4 = -1/30, f'''(n) = -n^-2. So coeff is ((-1/30)/24) * (-1) = 1/720.
# The term for n^-4 comes from B_6 / 6! * f^(5)(n).
# B_6 = 1/42, f^(5)(n) = -6n^-4. So coeff is ((1/42)/720) * (-6) = -1/5040.

# Define the numbers in the equation
c1 = 720
c2 = 5040

# Construct and print the formula for P(n)
# The formula represents the exponential of the correction terms.
formula_string = f"P(n) = exp(1/({c1}*n^2) - 1/({c2}*n^4))"

print(formula_string)
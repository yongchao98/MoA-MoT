import fractions

# The dimension-like exponent of the Frostman measure
s = fractions.Fraction(8, 5)

# According to a theorem by T. Wolff, the L^2 integral of the Fourier transform
# over a circle of radius r decays with an exponent of (1-s).
# The L^2 norm (which is the square root) decays with half of that exponent.
# So, the exponent c is given by the formula: c = (1 - s) / 2.
c = (1 - s) / 2

print(f"The Frostman measure exponent is s = {s.numerator}/{s.denominator} = {float(s)}")
intermediate_exponent = 1 - s
print(f"The exponent for the L^2 integral is 1 - s = 1 - {s} = {intermediate_exponent}")
print(f"The smallest possible c is (1 - s) / 2 = {intermediate_exponent} / 2 = {c}")
print(f"The value of c as a decimal is {float(c)}")

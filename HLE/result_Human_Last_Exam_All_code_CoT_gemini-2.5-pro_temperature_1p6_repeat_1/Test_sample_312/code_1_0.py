import fractions

# This script calculates the critical exponent 'c' for the decay of the Fourier
# transform of a Frostman measure, based on a well-known theorem in harmonic analysis.

# --- Problem Parameters ---
# The dimension of the space is R^d
d = 2
# The measure is an s-Frostman measure. s is given as a fraction.
s_numerator = 8
s_denominator = 5
# Divisor from the formula for the exponent
divisor = 2

# Create a fraction object for s for precise arithmetic
s = fractions.Fraction(s_numerator, s_denominator)

# --- Theoretical Background ---
print("The problem asks for the smallest possible value of c for the following relation:")
print("L^2 norm of \hat{\mu}(r\sigma) = O_{\epsilon}(r^{c+\epsilon})")
print("for any (8/5)-Frostman measure mu on R^2.\n")
print("According to a sharp theorem in Fourier analysis (by T. Wolff), the critical exponent 'c' is given by the formula:")
print("c = -(d - s) / 2")
print("where d is the dimension of the space and s is the Frostman exponent of the measure.\n")


# --- Calculation ---
print("We will now compute 'c' using the provided values.")

# Per the user request, we first output the numbers used in the final equation.
print("The numbers that constitute the equation c = -(d - s) / 2 are:")
print(f"d = {d}")
print(f"s = {s_numerator}/{s_denominator}")
print(f"The divisor in the formula is {divisor}\n")

# Check if the condition for the theorem (s > d/2) is met
if s > d / 2:
    print(f"The condition s > d/2 is satisfied, since s = {float(s):.2f} which is greater than d/2 = {d/2:.2f}.\n")
else:
    print(f"Note: The condition s > d/2 is not satisfied, which would make the theorem inapplicable.\n")


# Perform the step-by-step calculation
print("Calculation steps:")
d_minus_s = d - s
numerator_c = -d_minus_s
c = numerator_c / divisor

# Using the `fractions` module to show the arithmetic clearly
print(f"c = -({d} - {s}) / {divisor}")
print(f"  = -({fractions.Fraction(d)} - {s}) / {divisor}")
print(f"  = -({d_minus_s}) / {divisor}")
print(f"  = {numerator_c} / {divisor}")
print(f"  = {c}\n")

print(f"The smallest possible value for c is {c}, or {float(c)} as a decimal.")

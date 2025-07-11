from fractions import Fraction

# Step 1: Define the parameters given in the problem.
# The ambient space is R^n.
n = 2
# The measure is an alpha-Frostman measure.
alpha = Fraction(8, 5)

print(f"The problem is set in R^n with n = {n}.")
print(f"The measure is an alpha-Frostman measure with alpha = {alpha.numerator}/{alpha.denominator} = {float(alpha)}.")

# Step 2: The sharp decay exponent is known for alpha > 3/2. We check this condition.
condition = Fraction(3, 2)
print(f"\nA key theorem by Wolff gives the sharp decay exponent for alpha > 3/2.")
print(f"We check if our alpha satisfies this condition: {float(alpha)} > {float(condition)} is {alpha > condition}.")
if not alpha > condition:
    print("The condition is not met, the formula may not apply.")
else:
    print("The condition is met.")

# Step 3: Calculate the decay exponent for the L^2 integral.
# The decay exponent for the integral I(r) = \int_{S^1} |\hat{\mu}(r\sigma)|^2 d\sigma
# is given by the formula E = -(n - alpha/2). In R^2, this is E = -(2 - alpha/2).
print("\nThe exponent for the L^2 integral is 'E'. We calculate it as follows:")
E = -(n - alpha / 2)

print(f"E = -(n - alpha/2)")
print(f"E = -({n} - ({alpha.numerator}/{alpha.denominator}) / 2)")
print(f"E = -({n} - {alpha.numerator}/{alpha.denominator * 2})")
# To show the subtraction with a common denominator:
print(f"E = -({Fraction(n * 5, 5)} - {Fraction(4, 5)})")
print(f"E = -({Fraction(10 - 4, 5)})")
print(f"E = {E}")

# Step 4: Calculate the decay exponent for the L^2 norm.
# The L^2 norm is the square root of the integral, so its decay exponent 'c' is E/2.
print("\nThe exponent 'c' for the L^2 norm is half of 'E'.")
c = E / 2
print(f"c = E / 2")
print(f"c = ({E}) / 2")
print(f"c = {c}")

print("\n-------------------------------------------")
print(f"The smallest possible value for c is {c.numerator}/{c.denominator}, which is {float(c)}.")
print("-------------------------------------------")
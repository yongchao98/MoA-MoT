from fractions import Fraction

# Step 1: Define the function f(x) based on our derived form
# The roots are 0, 1, and alpha.
# f(x) = x * (x - 1) * (x - alpha)
# We found alpha = -5/8.
alpha = Fraction(-5, 8)

def f(x):
  # Convert x to a Fraction for exact calculations
  x = Fraction(x)
  return x * (x - 1) * (x + alpha)

# Step 2: Compute f(3)
x_val = 3
result = f(x_val)

# Step 3: Print the calculation step-by-step
# The final equation is f(3) = 3 * (3 - 1) * (3 + (-5/8))
term1 = Fraction(x_val)
term2 = Fraction(x_val - 1)
term3 = Fraction(x_val) + alpha

print("The function is f(x) = x(x-1)(x - 5/8).")
print("We need to compute f(3).")
print("f(3) = 3 * (3 - 1) * (3 - 5/8)")
print("The equation is composed of the following numbers:")
print(term1.numerator if term1.denominator == 1 else term1)
print(term2.numerator if term2.denominator == 1 else term2)
print(f"{term3.numerator}/{term3.denominator}")
print("The result of the multiplication is:")
print(f"{result.numerator}/{result.denominator}")

from fractions import Fraction

# Step 1: Define the roots based on the derivation.
# The roots are n1=0, n2=1, and r=-5/8.
r1 = Fraction(-5, 8)
r2 = 0
r3 = 1

# Step 2: Calculate the coefficients a, b, c from the roots using Vieta's formulas.
# f(x) = (x-r1)(x-r2)(x-r3) = x^3 - (r1+r2+r3)x^2 + (r1r2+r2r3+r3r1)x - r1r2r3
a = -(r1 + r2 + r3)
b = r1*r2 + r2*r3 + r3*r1
c = -(r1*r2*r3)

# Step 3: Define the function f(x)
def f(x_val):
    x = Fraction(x_val)
    return x**3 + a*x**2 + b*x + c

# Step 4: Compute the value of f(3)
x_to_eval = 3
result = f(x_to_eval)

# Step 5: Print the final equation and result
# The problem asks for the equation and the exact value.
# f(3) = 3^3 + a*3^2 + b*3 + c
term1 = 3**3
term2_val = a * 3**2
term3_val = b * 3
term4_val = c

print(f"The function is f(x) = x^3 + ({a})x^2 + ({b})x + ({c})")
print(f"We want to compute f(3).")
print(f"f(3) = 3^3 + ({a})*3^2 + ({b})*3 + ({c})")
print(f"f(3) = {term1} + ({term2_val}) + ({term3_val}) + ({term4_val})")
print(f"f(3) = {result}")

# The final answer in the required format
# <<<answer>>>
final_answer = result
# <<<87/4>>>
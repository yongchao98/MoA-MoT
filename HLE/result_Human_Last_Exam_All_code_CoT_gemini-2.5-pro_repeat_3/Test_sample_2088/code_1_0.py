import math

# The original expression is (12)^4 * (I_1 - I_2)^4, where I_1 and I_2 are the two integrals.
# After simplification and analytical integration, the value of the integral difference I = I_1 - I_2
# is found to be (5^(1/4)) / 12.

# Let's define the constants from the expression.
c1 = 12
power = 4

# Let's define the components of the integral's value.
# The value is integral_numerator / integral_denominator
integral_numerator_base = 5
integral_numerator_power_num = 1
integral_numerator_power_den = 4
integral_denominator = 12

# The overall expression is (c1)^power * (integral_value)^power, which can be simplified to (c1 * integral_value)^power.
# We will show the calculation step-by-step.

print(f"The expression to compute is ({c1})^({power}) * (Integral)^({power})")
print(f"The value of the integral difference is ({integral_numerator_base}^({integral_numerator_power_num}/{integral_numerator_power_den})) / {integral_denominator}")
print(f"So we compute: ({c1})^({power}) * (({integral_numerator_base}^({integral_numerator_power_num}/{integral_numerator_power_den})) / {integral_denominator})^({power})")
print(f"= ({c1} * (({integral_numerator_base}^({integral_numerator_power_num}/{integral_numerator_power_den})) / {integral_denominator}))^{power}")
print(f"= (({c1}/{integral_denominator}) * {integral_numerator_base}^({integral_numerator_power_num}/{integral_numerator_power_den}))^{power}")
print(f"= (1 * {integral_numerator_base}^({integral_numerator_power_num}/{integral_numerator_power_den}))^{power}")
print(f"= ({integral_numerator_base}^({integral_numerator_power_num}/{integral_numerator_power_den}))^{power}")
print(f"= {integral_numerator_base}^(({integral_numerator_power_num}/{integral_numerator_power_den}) * {power})")

final_power = (integral_numerator_power_num / integral_numerator_power_den) * power
print(f"= {integral_numerator_base}^{int(final_power)}")

final_result = integral_numerator_base ** final_power
print(f"= {int(final_result)}")
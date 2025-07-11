import math

# The value of the integral I is (5^(1/4))/12
c1 = 12
c2 = 5
integral_val_numerator = c2**(1/4)
integral_val_denominator = c1
integral_val = integral_val_numerator / integral_val_denominator

# The original expression is (12**4) * (integral_val**4)
base_val = 12
exponent = 4

# Let's compute the result step-by-step and print the equation
term1 = base_val**exponent
term2_num = integral_val_numerator**exponent
term2_den = integral_val_denominator**exponent

final_result = term1 * (term2_num / term2_den)

print(f"The integral evaluates to I = ( {c2}^(1/4) ) / {c1}")
print(f"The expression to compute is ( {base_val}^{exponent} ) * I^{exponent}")
print(f"This is equal to ( {term1} ) * ( ( {c2}^(1/4) )^{exponent} / {c1}^{exponent} )")
print(f"Which simplifies to ( {term1} ) * ( {int(term2_num)} / {term2_den} )")
print(f"The final result is {int(final_result)}")

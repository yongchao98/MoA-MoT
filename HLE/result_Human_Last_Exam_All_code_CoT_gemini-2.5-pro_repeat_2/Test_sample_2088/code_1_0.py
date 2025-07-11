import math

# Define the components of the final calculation
c1 = 12
p1 = 4
integral_val_num_base = 5
integral_val_num_power_numerator = 1
integral_val_num_power_denominator = 4
integral_val_den = 12

# Perform the final calculation
# Expression: (c1^p1) * ( (integral_val_num_base^(1/4)) / integral_val_den )^p1
result = (c1**p1) * ((integral_val_num_base**(1/4)) / integral_val_den)**p1

print("After simplifying the integral expression, its value is 5^(1/4) / 12.")
print("We substitute this back into the original expression to get the final equation:")
print(f"({c1})^{p1} * (({integral_val_num_base}^({integral_val_num_power_numerator}/{integral_val_num_power_denominator})) / {integral_val_den})^{p1}")
print("\nThis simplifies as follows:")
print(f"= ({c1})^{p1} * (({integral_val_num_base}^({integral_val_num_power_numerator * p1}/{integral_val_num_power_denominator})) / ({integral_val_den})^{p1})")
print(f"= ({c1})^{p1} * ({integral_val_num_base}^1 / {integral_val_den}^{p1})")
print(f"= (({c1})^{p1} / ({integral_val_den})^{p1}) * {integral_val_num_base}")
print(f"= 1 * {integral_val_num_base}")
print(f"= {int(result)}")

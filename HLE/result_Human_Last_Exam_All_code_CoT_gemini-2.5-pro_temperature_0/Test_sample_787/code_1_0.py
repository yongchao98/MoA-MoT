# This script calculates the limit of the sequence g_n based on the derivation above.

# The limit L is of the form 2^a * 3^b * 5^c.
# From the mathematical analysis, we found the exponents to be:
a = 10  # v_2(L)
b = 2   # v_3(L)
c = 1   # v_5(L)

# Calculate the components of the final product
power_of_2 = 2**a
power_of_3 = 3**b
power_of_5 = 5**c

# Calculate the final limit
limit_value = power_of_2 * power_of_3 * power_of_5

# Print the components of the final equation as requested
print(f"The limit is calculated as 2^{a} * 3^{b} * 5^{c}")
print(f"2^{a} = {power_of_2}")
print(f"3^{b} = {power_of_3}")
print(f"5^{c} = {power_of_5}")
print(f"The final result is {power_of_2} * {power_of_3} * {power_of_5} = {limit_value}")

print(f"\n{limit_value}")
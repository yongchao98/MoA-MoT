import math

# The dominant root of the characteristic equation x^4 - 88x - 57 = 0 is lambda = 2 + sqrt(7).
# The limit is L = 2 * ln(lambda).
# We need to compute the integer part of 10^4 * L.

dom_root_lambda = 2 + math.sqrt(7)
limit_L = 2 * math.log(dom_root_lambda)
value = 10000 * limit_L

# The final result is the integer part of the calculated value.
integer_part = int(value)

print(f"The characteristic equation is x^4 - 88x - 57 = 0.")
print(f"The roots are 2+sqrt(7), 2-sqrt(7), -2+i*sqrt(15), -2-i*sqrt(15).")
print(f"The dominant root is lambda = 2 + sqrt(7) = {dom_root_lambda}")
print(f"The limit is L = 2 * ln(lambda) = {limit_L}")
print(f"The value to compute is 10^4 * L = {value}")
print(f"The integer part of the value is {integer_part}.")
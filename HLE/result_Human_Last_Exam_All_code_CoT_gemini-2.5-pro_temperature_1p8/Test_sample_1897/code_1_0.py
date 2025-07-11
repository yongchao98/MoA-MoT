import math

# The derived closed-form expression for the infinite product is 2**(1/8) * e**(-pi/24).
# This script computes its numerical value and displays the components of the expression.

# The numbers in the final expression: 2**(1/8) * e**(-pi/24)
base1 = 2
exponent1_num = 1
exponent1_den = 8

base2 = math.e
exponent2_num = -1
exponent2_den = 24
pi = math.pi

# Calculate the numerical value
value = (base1**(exponent1_num / exponent1_den)) * (base2**((exponent2_num * pi) / exponent2_den))

# Print the components and the final result
print(f"The closed-form expression is: {base1}**({exponent1_num}/{exponent1_den}) * e**({exponent2_num}*π/{exponent2_den})")
print(f"The numerical value of the product is: {value}")
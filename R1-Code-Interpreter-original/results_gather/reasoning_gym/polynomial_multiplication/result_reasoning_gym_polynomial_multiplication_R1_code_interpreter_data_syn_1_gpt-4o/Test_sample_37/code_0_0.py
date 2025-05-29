# Define the coefficients
a = 22
b = -42
c = 31
d = -51

# Calculate each part of the expression
first = a * c
outer = a * d
inner = b * c
last = b * d

# Combine the results
simplified_expression = f"{first}*z + ({outer} + {inner})*z + {last}"
print(simplified_expression)
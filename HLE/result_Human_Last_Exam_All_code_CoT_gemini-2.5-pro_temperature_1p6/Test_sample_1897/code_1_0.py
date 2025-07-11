import math

# The closed-form expression is 2^(1/8) * e^(-pi/24)
base1 = 2
exponent1 = 1/8
base2 = math.e
exponent2 = -math.pi/24

# Calculate the final value
result = (base1 ** exponent1) * (base2 ** exponent2)

print("The closed expression for the infinite product is: A ** B * C ** D")
print("Where the numbers in the equation are:")
print(f"A = {base1}")
print(f"B = {exponent1}")
print(f"C = {base2}")
print(f"D = {exponent2}")
print("-" * 20)
print(f"The final numerical value is: {result}")

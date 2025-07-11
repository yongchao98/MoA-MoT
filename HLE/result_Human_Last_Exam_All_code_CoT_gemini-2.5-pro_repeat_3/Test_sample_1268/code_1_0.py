import math

# The problem asks for the relationship between the upper bound for the
# norm (k_k,inf) and the covolume (V) for number fields related to
# squarefree integers. We derive this for imaginary quadratic fields.
# The relationship is k_k,inf <= (4/pi) * V.

# Define the numbers in the equation's constant.
numerator = 4
denominator_symbol = "pi"
denominator_value = math.pi

# Calculate the value of the constant.
constant_value = numerator / denominator_value

print("The relationship is of the form: k_k,inf <= C * V")
print("The constant C is derived from the following numbers:")
print(f"Numerator: {numerator}")
print(f"Denominator: {denominator_symbol} (value ≈ {denominator_value})")
print(f"\nThe resulting constant C is {numerator}/{denominator_symbol} ≈ {constant_value}")

print("\nTherefore, the final equation showing each number is:")
# The f-string below prints the full equation with its components.
print(f"k_k,inf <= ({numerator} / {denominator_value}) * V")
print("\nWhich simplifies to:")
print(f"k_k,inf <= {constant_value} * V")
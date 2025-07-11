import math

# The analytical solution to the integral is (108 * pi) / 35.
# Here are the components of that final equation:
numerator_coefficient = 108
denominator = 35
pi_value = math.pi

# Calculate the final numerical result
final_value = (numerator_coefficient * pi_value) / denominator

print("The integral of f(x,y)=z^2(x^2+y^2) inside the cone is given by the equation (108 * pi) / 35.")
print("\nHere are the numbers in the final equation:")
print(f"Numerator Coefficient: {numerator_coefficient}")
print(f"Pi: {pi_value}")
print(f"Denominator: {denominator}")
print(f"\nThe final equation with the computed value is:")
print(f"({numerator_coefficient} * {pi_value}) / {denominator} = {final_value}")
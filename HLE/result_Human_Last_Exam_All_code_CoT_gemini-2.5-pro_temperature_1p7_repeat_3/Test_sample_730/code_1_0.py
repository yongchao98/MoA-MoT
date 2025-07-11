# Define the integer components for the calculation
density_numerator = 9
density_denominator = 10
radius_numerator = 1
radius_denominator = 2
pi_approximation = 3
volume_constant_numerator = 4
volume_constant_denominator = 3

# Calculate the numerator and denominator of the result
result_numerator = density_numerator * volume_constant_numerator * pi_approximation * (radius_numerator**3)
result_denominator = density_denominator * volume_constant_denominator * (radius_denominator**3)

# To simplify the final fraction, we find the greatest common divisor
def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

common_divisor = gcd(result_numerator, result_denominator)
simplified_numerator = result_numerator // common_divisor
simplified_denominator = result_denominator // common_divisor

# Print the final equation, showing each number used in the calculation explicitly, as requested.
print(
    f"mass = {density_numerator}/{density_denominator} "
    f"* {volume_constant_numerator}/{volume_constant_denominator} "
    f"* {pi_approximation} "
    f"* {radius_numerator}/{radius_denominator} "
    f"* {radius_numerator}/{radius_denominator} "
    f"* {radius_numerator}/{radius_denominator} "
    f"= {simplified_numerator}/{simplified_denominator}"
)
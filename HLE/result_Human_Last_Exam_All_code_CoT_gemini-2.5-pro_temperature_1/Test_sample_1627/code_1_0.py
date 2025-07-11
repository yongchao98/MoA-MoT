import sympy

# Define the symbols
# Pi: Internal Pressure
# a: Inner radius
# b: Outer radius
Pi, a, b = sympy.symbols('Pi a b')

# Lamé's equation for maximum tangential stress (hoop stress) at the inner radius (r=a)
# under internal pressure Pi and zero external pressure.
# sigma_t_max = Pi * (b**2 + a**2) / (b**2 - a**2)
# We want to find the ratio: sigma_t_max / Pi
ratio_expr = (b**2 + a**2) / (b**2 - a**2)

print("The formula for the ratio of maximum tangential stress to internal pressure is:")
print(f"Ratio = (b^2 + a^2) / (b^2 - a^2)")
print("\nGiven that the outer radius (b) is twice the inner radius (a), we have b = 2a.")

# Substitute b = 2a into the expression
ratio_with_substitution = ratio_expr.subs(b, 2*a)

print("\nSubstituting b = 2a into the formula:")
# To show the substitution clearly, we represent the terms as strings
numerator_str = f"((2a)^2 + a^2)"
denominator_str = f"((2a)^2 - a^2)"
print(f"Ratio = {numerator_str} / {denominator_str}")

# Simplify the expression
simplified_numerator_str = f"(4a^2 + a^2)"
simplified_denominator_str = f"(4a^2 - a^2)"
print(f"Ratio = {simplified_numerator_str} / {simplified_denominator_str}")

final_numerator_str = "5a^2"
final_denominator_str = "3a^2"
print(f"Ratio = {final_numerator_str} / {final_denominator_str}")

# Final simplification by canceling a^2
final_ratio = sympy.simplify(ratio_with_substitution)
numerator, denominator = final_ratio.as_numer_denom()

print(f"\nAfter canceling a^2, the final ratio is:")
print(f"Ratio = {numerator}/{denominator}")

# Also print as a decimal
print(f"As a decimal, the ratio is approximately: {float(final_ratio):.4f}")
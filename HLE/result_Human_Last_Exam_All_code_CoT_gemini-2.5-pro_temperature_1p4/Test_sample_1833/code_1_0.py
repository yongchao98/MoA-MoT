import math

# Define the given variables
m = 1.0  # kg, mass of the ring
M = 1.0  # kg, mass of the object
g = 9.8  # m/s^2, acceleration due to gravity
theta_deg = 60.0  # degrees

# Convert angle to radians for trigonometric functions
theta_rad = math.radians(theta_deg)
sin_theta = math.sin(theta_rad)

# For the case m=M=1 and theta=60 degrees, the derived formula simplifies significantly.
# T = (94/25) * M * g * sin(theta)
# Let's show the calculation with the specific numbers.
factor_numerator = 94
factor_denominator = 25

# Calculate the tension
tension = (factor_numerator / factor_denominator) * M * g * sin_theta

# Output the equation and the final result
print("The simplified equation for tension in this specific case (m=M, theta=60 deg) is:")
print(f"T = ({factor_numerator}/{factor_denominator}) * M * g * sin(theta)")
print("\nSubstituting the given values:")
print(f"T = ({factor_numerator}/{factor_denominator}) * {M} * {g} * sin({theta_deg} deg)")
print(f"T = {factor_numerator/factor_denominator} * {M} * {g} * {sin_theta:.3f}")
print(f"\nThe calculated tension is {tension:.2f} Newtons.")

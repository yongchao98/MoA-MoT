import math

# This script calculates the maximum achievable ratio of bidirectional conical power
# to line intensity for the described rotating sphere.

# Plan steps implemented in code:
# 1. The angular dependence of the magnetic dipole field intensity is defined as
#    I(θ) = 3*cos²(θ) + 1.
# 2. We find the minimum value of this intensity (the denominator).
# 3. We calculate the maximum value of the conical power (the numerator) by
#    integrating I(θ) over the specified bidirectional cone.
# 4. We compute and display the final ratio, outputting each number involved.

# --- Calculation ---

# DENOMINATOR: Minimum Line Intensity
# The intensity I(θ) is minimized when cos(θ) = 0, which occurs at θ = π/2.
theta_for_min_intensity = math.pi / 2
min_line_intensity = 3 * (math.cos(theta_for_min_intensity)**2) + 1

# NUMERATOR: Maximum Conical Power
# P_max = ∫[0 to 2π]dφ ∫[cone]dθ sin(θ) * (3*cos²(θ) + 1)
# The integral over φ gives a factor of 2π.
# For the θ integral, we find the anti-derivative of (3*cos²(θ) + 1)*sin(θ).
# Using substitution u = cos(θ), we integrate -(3u² + 1)du, which gives -u³ - u.
# So, the anti-derivative F(θ) is -cos³(θ) - cos(θ).

def anti_derivative(theta):
    """The anti-derivative F(θ) of the integrand."""
    cos_theta = math.cos(theta)
    return - (cos_theta**3) - cos_theta

# The definite integral over the first cone (θ from 0 to π/4) is F(π/4) - F(0).
integral_first_cone = anti_derivative(math.pi / 4) - anti_derivative(0)

# The definite integral over the second cone (θ from 3π/4 to π) is F(π) - F(3π/4).
integral_second_cone = anti_derivative(math.pi) - anti_derivative(3 * math.pi / 4)

# The total conical power is 2π times the sum of the results of the definite integrals.
max_conical_power = 2 * math.pi * (integral_first_cone + integral_second_cone)

# The final maximum ratio is the power divided by the intensity.
max_ratio = max_conical_power / min_line_intensity

# --- Outputting the Final Equation ---

print("The maximum achievable ratio is determined by the equation: Ratio = P_max / I_min")
print("\nThe values for each number in this final equation are:")
print(f"P_max (Maximum Bidirectional Conical Power) = {max_conical_power}")
print(f"I_min (Minimum Line Intensity)             = {min_line_intensity}")
print(f"Maximum Ratio = {max_conical_power} / {min_line_intensity} = {max_ratio}")
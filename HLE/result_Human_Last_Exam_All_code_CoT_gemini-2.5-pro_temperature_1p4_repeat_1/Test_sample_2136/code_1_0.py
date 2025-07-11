# Turbulent plasma velocity fluctuations problem

# The problem asks for the value of the integral ∫(∂u/∂t)² dx.
# For the given PDE, there is a remarkable identity that relates this integral
# to the spatial gradient at a stationary point of the solution.
# The identity is: I = (144/5) * |∂u/∂x|^5
# We are given that ∂u/∂x at the stationary point (x_0, τ) is -1.

# Coefficients and given values
numerator = 144
denominator = 5
gradient = -1
power = 5

# Calculate the value of the integral
# We use abs() for the absolute value, as is proper for the identity.
result = (numerator / denominator) * (abs(gradient) ** power)

# Print the final equation with all numbers
print(f"The calculation is based on the identity: ∫(∂u/∂t)² dx = ({numerator}/{denominator}) * |∂u/∂x(x₀,τ)|⁵")
print(f"Given ∂u/∂x(x₀,τ) = {gradient}, the equation becomes:")
print(f"Result = ({numerator}/{denominator}) * abs({gradient}) ** {power}")
print(f"Result = {result}")

# The separatrix is a curve in the (u, d) phase plane that separates
# trajectories with different behaviors. Based on the analysis of the system's
# dynamics, the separatrix is found to be a parabola.

# The equation of the separatrix is d = -u^2.
# To present this in the format d + f(u) = 0, we can write it as d + u^2 = 0.
# The following code prints this equation, explicitly stating the coefficients and powers.

# Coefficients and powers in the equation: 1*d + 1*u^2 = 0
coeff_d = 1
coeff_u_sq = 1
power_u = 2
constant = 0

# Print the equation with all its numbers
print("The equation for the separatrix is:")
print(f"{coeff_d}*d + {coeff_u_sq}*u**{power_u} = {constant}")

import math

# The problem asks for an expression for the angle θ at which the rod begins to slide.
# Based on the physics derivation, this angle depends only on the coefficient of friction μ.

# The final relationship is tan(θ) = μ.
# Or, expressing θ explicitly: θ = arctan(μ).

# The following code prints this final expression.
# We use unicode characters for theta and mu for better readability.

theta_symbol = "θ"
mu_symbol = "μ"

# Print the final equation with each component part by part
print("The final expression for the angle is:")
print(theta_symbol, end="")
print(" = ", end="")
print(f"arctan({mu_symbol})")

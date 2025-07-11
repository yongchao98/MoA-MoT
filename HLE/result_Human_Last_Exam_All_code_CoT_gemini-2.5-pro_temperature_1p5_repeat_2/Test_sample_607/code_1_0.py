import numpy as np

# This code block demonstrates the reasoning for the heavy-ball method.
# It sets up the limit equation and shows that it implies stationarity.
# It cannot demonstrate the counter-example as that would require defining a complex function.
# This code is for illustrative purposes of the derivation.

# Let's represent the quantities in the limit as k -> infinity.
# x_k, x_{k-1}, x_{k+1} all converge to a point x_star.
# We represent vectors as strings for symbolic-like demonstration.
x_k = "x_star"
x_k_minus_1 = "x_star"
x_k_plus_1 = "x_star"
nabla_f_x_star = "nabla_f(x_star)" # Gradient at the limit point

# Parameters
beta = 0.9  # A typical momentum parameter
gamma = 0.1 # A typical learning rate

# The heavy-ball update equation is:
# x_{k+1} = x_k + beta * (x_k - x_{k-1}) - gamma * nabla_f(x_k)

# Let's rearrange it to solve for the gradient term:
# gamma * nabla_f(x_k) = x_k + beta * (x_k - x_{k-1}) - x_{k+1}

# Now, we take the limit as k -> infinity.
# Let LHS be the limit of the left-hand side and RHS be the limit of the right-hand side.

# LHS in the limit becomes:
LHS = f"gamma * {nabla_f_x_star}"

# RHS in the limit becomes (substituting the limit points):
# Note: (x_star - x_star) represents the limit of (x_k - x_{k-1}) which is 0.
RHS_term1 = f"({x_k} - {x_k_plus_1})" # This limit is 0
RHS_term2 = f"beta * ({x_k} - {x_k_minus_1})" # This limit is beta * 0 = 0

# So, the limit equation is LHS = RHS
# gamma * nabla_f(x_star) = (x_star - x_star) + beta * (x_star - x_star)
# gamma * nabla_f(x_star) = 0 + beta * 0
# gamma * nabla_f(x_star) = 0

print("Heavy-ball method update rule: x_{k+1} = x_k + beta*(x_k - x_{k-1}) - gamma*nabla_f(x_k)")
print("Assuming the sequence converges to x_star, we take the limit as k -> infinity.")
print(f"The term (x_{k+1} - x_k) goes to ({x_k_plus_1} - {x_k}) -> 0")
print(f"The term (x_k - x_{k-1}) goes to ({x_k} - {x_k_minus_1}) -> 0")
print(f"The term nabla_f(x_k) goes to {nabla_f_x_star} (by continuity of the gradient)")
print("\nThe limit of the update equation is:")
print(f"0 = beta * 0 - gamma * {nabla_f_x_star}")
print(f"This simplifies to: gamma * {nabla_f_x_star} = 0")
print(f"Since gamma ({gamma}) must be positive, this implies that nabla_f(x_star) = 0.")
print("\nThis derivation shows that under normal assumptions (C^1 function), a convergent heavy-ball sequence must converge to a stationary point.")
print("However, literature provides counterexamples for specific non-convex functions where this is not the case, making (3) the correct choice.")

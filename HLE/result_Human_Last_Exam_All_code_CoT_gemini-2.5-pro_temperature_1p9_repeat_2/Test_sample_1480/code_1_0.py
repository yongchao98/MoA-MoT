# Given values
phi = 0.7
mu = 0.9

# The rule activation level is calculated using a t-norm,
# typically the minimum operator.
activation_level = min(phi, mu)

# Print the equation and the result
print(f"Rule Activation Level = t-norm(Phi_k(x'), μG_k_j(y_j))")
print(f"Rule Activation Level = min({phi}, {mu})")
print(f"Rule Activation Level = {activation_level}")

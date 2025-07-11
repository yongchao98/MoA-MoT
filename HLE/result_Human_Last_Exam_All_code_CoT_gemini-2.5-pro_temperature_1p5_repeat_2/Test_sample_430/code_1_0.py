import numpy as np
import sys

# Step 1 & 2: Define parameters based on the problem description and the k=1 assumption.
# The parameters are given in dimensionless units of k_B*T, so beta = 1/(k_B*T) is already incorporated.
# For example, beta_eps1 represents beta * eps1.
# We assume k=1, so eps_l / (k_B*T) = (0.02)^1 = 0.02.
beta_eps1 = 0.1
beta_eps_l = 0.02
beta_mu = 0.15
z_l = 4

# The coordination number for inter-layer interaction, z_inter, is not used in the k=1 model.
# The temperature T=318K is also not needed as all energies are given relative to k_B*T.

# Step 3: Define the function for the self-consistency equation.
# The equation is: theta = 1 / (1 + exp(beta * (E_eff - mu)))
# where E_eff = eps1 + z_l * eps_l * theta.
# Let's define the right-hand side as a function f(theta).
def f(theta, beta_mu_val, beta_eps1_val, beta_eps_l_val, z_l_val):
    """
    Calculates the right-hand side of the self-consistency equation.
    theta = 1 / (1 + exp(beta*eps1 + beta*z_l*eps_l*theta - beta*mu))
    """
    exponent = beta_eps1_val + z_l_val * beta_eps_l_val * theta - beta_mu_val
    return 1.0 / (1.0 + np.exp(exponent))

# Step 4: Solve the equation theta = f(theta) using fixed-point iteration.
theta = 0.5  # Initial guess for the coverage
# Iterate up to 100 times to find the solution. Convergence is usually fast.
for _ in range(100):
    theta_new = f(theta, beta_mu, beta_eps1, beta_eps_l, z_l)
    # Check for convergence
    if np.abs(theta_new - theta) < 1e-9:
        theta = theta_new
        break
    theta = theta_new

# Step 5: Output the explanation and the final result.
print("Derivation and Method")
print("=====================")
print("The problem is ambiguous as some parameters ('k', 'eps_inter') are not specified. We proceed by assuming the problem is intended to be solved for the case of monolayer adsorption (k=1), as this provides a self-consistent problem that utilizes most of the given parameters.")
print("\nFor a monolayer (k=1), the average number of layers per site is the fractional coverage, θ.")
print("Using the mean-field approximation, the coverage θ must satisfy the self-consistency equation:")
print("\n  θ = 1 / (1 + exp(β * (E_eff - μ)))")
print("\nwhere β=1/(k_B*T) and the effective energy E_eff is:")
print("\n  E_eff = ε₁ + z_ℓ * ε_ℓ * θ")
print("\nSubstituting the given numerical values (already scaled by β), the equation to solve is:")
# This part fulfills the requirement to output each number in the final equation.
print(f"\n  θ = 1 / (1 + exp( ({beta_eps1}) + ({z_l}) * ({beta_eps_l}) * θ - ({beta_mu}) ))")

print("\n\nResult")
print("========")
print("This transcendental equation is solved numerically via fixed-point iteration.")
print("The calculated average number of adsorbed layers per site is:")
# sys.stdout.write ensures only the final number is available for capture.
sys.stdout.write(f"<<<{theta:.4f}>>>\n")

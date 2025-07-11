import math

# Step 1 & 2: Define parameters
# The physical model is a lattice gas under the mean-field approximation.
# The total coordination number z is the sum of the horizontal and vertical parts.
z_horizontial = 4
z_vertical = 8
z = z_horizontial + z_vertical

# Step 3 & 4: Define dimensionless energy parameters from the problem statement.
# The temperature T = 300 K is implicit in these ratios.
# beta * mu = mu / (k_B * T)
beta_mu = 0.1
# beta * epsilon = epsilon / (k_B * T)
beta_epsilon = -1 / (2 * math.pi)

# The self-consistency equation is:
# <n> = 1 / (exp(z * beta_epsilon * <n> - beta_mu) + 1)
# Let's define the coefficients for the exponent term.
coeff_n = z * beta_epsilon
constant_term = -beta_mu

# Print the numerical equation being solved, as requested.
print("The self-consistency equation for the average occupancy <n> is:")
print(f"⟨n⟩ = 1 / (exp({z} * ({beta_epsilon:.4f}) * ⟨n⟩ - {beta_mu:.1f}) + 1)")
print(f"Which simplifies to:")
print(f"⟨n⟩ = 1 / (exp({coeff_n:.4f} * ⟨n⟩ - {beta_mu:.1f}) + 1)")
print("-" * 30)

# Step 5: Solve the equation for <n> using an iterative method.
# We start with an initial guess for n_avg (<n>).
n_avg = 0.5
# We iterate a fixed number of times, which is sufficient for convergence.
for _ in range(100):
    exponent = coeff_n * n_avg + constant_term
    n_avg = 1 / (math.exp(exponent) + 1)

# Step 6: Print the final result.
print(f"The final calculated average occupancy per site is:")
print(f"⟨n⟩ = {n_avg:.3f}")

print(f"\n<<<{n_avg:.3f}>>>")
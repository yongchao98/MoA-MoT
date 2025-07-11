import numpy as np

# --- 1. Define Parameters ---
# The problem gives parameters in terms of kB*T, so we can work with
# dimensionless quantities beta*epsilon and beta*mu where beta = 1/(kB*T).

# Dimensionless interaction energy: beta * epsilon
beta_epsilon = -1.0 / (2.0 * np.pi)

# Dimensionless chemical potential: beta * mu
beta_mu = 0.1

# Coordination numbers
z_horizontial = 4
z_vertical = 8
# The total coordination number z is the sum of horizontal and vertical neighbors
# for a site not in the first layer or at the surface of the multilayer film.
z = z_horizontial + z_vertical

# --- 2. Solve the Self-Consistent Equation ---
# The mean-field theory leads to the self-consistent equation for the
# average occupancy <n>:
# <n> = 1 / (exp(beta * (z * epsilon * <n> - mu)) + 1)
# We solve this using a simple fixed-point iteration method: n_new = f(n_old).

n = 0.5  # Initial guess for the average occupancy <n>
# Iterate 100 times to ensure convergence to the stable solution
for _ in range(100):
    # This is the exponent in the Fermi-Dirac like distribution
    exponent = z * beta_epsilon * n - beta_mu
    n = 1.0 / (np.exp(exponent) + 1.0)

# --- 3. Output the Results ---
# The problem requires printing the final equation with numerical values inserted.

print("The self-consistent mean-field equation for the average occupancy <n> is:")
print("<n> = 1 / (exp(beta * (z * epsilon * <n> - mu)) + 1)")
print("\nPlugging in the dimensionless parameters:")

z_str = f"z = {z_horizontial} + {z_vertical} = {z}"
beta_epsilon_str = f"beta * epsilon = -1 / (2 * pi) = {beta_epsilon:.3f}"
beta_mu_str = f"beta * mu = {beta_mu}"

print(z_str)
print(beta_epsilon_str)
print(beta_mu_str)

# Constructing the final equation string with numerical values
final_equation_str = f"<n> = 1 / (exp({z} * ({beta_epsilon:.3f}) * <n> - {beta_mu}) + 1)"
print("\nThe specific equation to solve is:")
print(final_equation_str)

print(f"\nSolving this equation numerically gives the average occupancy per site:")
# Round the final answer to three decimal places
final_n = round(n, 3)
print(f"<n> = {final_n}")
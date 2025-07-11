import math

# Step 1: Define the parameters based on the problem statement.
# The chemical potential is given in units of k_B*T.
mu_over_kBT = 0.1

# The interaction energy is given in units of k_B*T.
epsilon_over_kBT = -1 / (2 * math.pi)

# The coordination numbers for the lattice.
z_horizontal = 4
z_vertical = 8

# In the mean-field approximation, we consider the total number of nearest neighbors.
z = z_horizontal + z_vertical

# Step 2: Formulate the self-consistency equation.
# The mean-field equation for the average occupancy <n> is:
# <n> = 1 / (exp((z * epsilon * <n> - mu) / (k_B * T)) + 1)
# We can rewrite this using the dimensionless parameters defined above:
# <n> = 1 / (exp(z * (epsilon/k_B*T) * <n> - (mu/k_B*T)) + 1)

print("This script solves the self-consistency equation for the average occupancy <n> in a mean-field lattice gas model.")
print("The equation to solve is: <n> = 1 / (exp(z * (epsilon/kBT) * <n> - (mu/kBT)) + 1)\n")
print("Substituting the given parameter values:")
print(f"z = {z}")
print(f"epsilon/kBT = {epsilon_over_kBT:.4f}")
print(f"mu/kBT = {mu_over_kBT:.4f}\n")
print("The final equation being solved is:")
print(f"<n> = 1 / (exp({z} * ({epsilon_over_kBT:.4f}) * <n> - {mu_over_kBT:.4f}) + 1)")
print(f"<n> = 1 / (exp({z*epsilon_over_kBT:.4f} * <n> - {mu_over_kBT:.4f}) + 1)\n")

# Step 3: Solve the equation numerically using an iterative method.
# We start with an initial guess for <n> (e.g., 0.5) and iterate until convergence.
n = 0.5  # Initial guess for the average occupancy
tolerance = 1e-8
max_iterations = 100

for i in range(max_iterations):
    # Calculate the exponent term
    exponent = z * epsilon_over_kBT * n - mu_over_kBT
    
    # Calculate the next value of n
    n_next = 1 / (math.exp(exponent) + 1)
    
    # Check for convergence
    if abs(n_next - n) < tolerance:
        break
    
    n = n_next

# Step 4: Print the final result.
final_occupancy = n
print(f"After solving numerically, the average occupancy per site is:")
print(f"<n> = {final_occupancy:.3f}")

# Final answer in the required format
print(f"<<<{final_occupancy:.3f}>>>")
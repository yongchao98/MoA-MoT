import math

# 1. Define the given dimensionless parameters from the problem.
# The temperature T and Boltzmann constant kB are not needed individually,
# as the relevant quantities are given in terms of kB*T.
mu_over_kBT = 0.1
epsilon_over_kBT = -1.0 / (2.0 * math.pi)
z_horizontal = 4
z_vertical = 8
z = z_horizontal + z_vertical

# 2. The self-consistency equation for the average occupancy <n> is:
#    ln(<n> / (1 - <n>)) = mu/(kB*T) - z*epsilon/(kB*T) * <n>
#    We can write this as: ln(<n> / (1 - <n>)) = C1 + C2 * <n>
C1 = mu_over_kBT
C2 = -z * epsilon_over_kBT

# 3. Print the final equation with its numerical coefficients as requested.
print("The self-consistency equation for the average occupancy <n> is:")
print(f"ln(<n> / (1 - <n>)) = {C1:.1f} - ({z}) * ({epsilon_over_kBT:.4f}) * <n>")
print(f"Which simplifies to:")
print(f"ln(<n> / (1 - <n>)) = {C1:.1f} + {C2:.4f} * <n>")

# 4. Solve the equation numerically using an iterative approach.
#    The iterative form is: <n>_{i+1} = 1 / (1 + exp(-(C1 + C2 * <n>_i)))
#    Initialize <n> with a starting guess.
n_avg = 0.5

# Iterate a sufficient number of times for the value to converge.
for _ in range(100):
    exponent = C1 + C2 * n_avg
    n_avg = 1.0 / (1.0 + math.exp(-exponent))

# 5. Print the final calculated value, rounded to three decimal places.
print("\nSolving this equation numerically yields:")
print(f"The average occupancy per site <n> is: {n_avg:.3f}")
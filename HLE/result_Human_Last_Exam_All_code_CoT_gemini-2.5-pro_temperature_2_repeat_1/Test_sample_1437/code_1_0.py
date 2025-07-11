import numpy as np

# This script calculates the dimensionless value of the Lindhard polarization function
# at zero frequency and momentum transfer for a 3D homogeneous electron gas.

# Plan:
# 1. Define physical constants. Since we are calculating a dimensionless quantity,
#    their specific values do not matter, and they will cancel out. We set them to 1.
# 2. Choose an arbitrary electron density `n`.
# 3. Calculate the Fermi wavevector `k_F` from `n`.
# 4. Calculate the Fermi energy `epsilon_F` from `k_F`.
# 5. Calculate the density of states at the Fermi energy `D_epsilon_F`.
# 6. The Lindhard function at this limit is `Pi_00 = -D_epsilon_F`.
# 7. Compute the dimensionless ratio `result = Pi_00 * epsilon_F / n`.

# Step 1: Define physical constants (we can set them to 1)
hbar = 1.0  # Reduced Planck constant
m = 1.0     # Electron mass

# Step 2: Choose an arbitrary electron density
n = 1.0     # Electron density (e.g., in electrons per cubic meter)

# Step 3: Calculate the Fermi wavevector k_F from the electron density n
# The formula for a 3D electron gas is n = k_F^3 / (3 * pi^2)
k_F = (3 * np.pi**2 * n)**(1/3)

# Step 4: Calculate the Fermi energy epsilon_F
# The formula is epsilon_F = (hbar^2 * k_F^2) / (2 * m)
epsilon_F = (hbar**2 * k_F**2) / (2 * m)

# Step 5: Calculate the total density of states at the Fermi energy, D(epsilon_F)
# The formula (including spin) is D(epsilon_F) = m * k_F / (pi^2 * hbar^2)
D_epsilon_F = (m * k_F) / (np.pi**2 * hbar**2)

# Step 6: The Lindhard polarization function at k=0, omega=0 is Pi(0,0) = -D(epsilon_F)
Pi_00 = -D_epsilon_F

# Step 7: To obtain a universal numerical constant, we compute the dimensionless quantity:
# result = Pi(0,0) / (n / epsilon_F) which is equivalent to Pi(0,0) * epsilon_F / n
result = Pi_00 * epsilon_F / n

# The prompt requires outputting the numbers in the final equation.
print("The final calculation is for the dimensionless quantity Pi(0,0) * epsilon_F / n.")
print(f"The Lindhard function value Pi(0,0) is: {Pi_00}")
print(f"The Fermi energy epsilon_F is: {epsilon_F}")
print(f"The electron density n is: {n}")
print(f"The full equation with these values is: ({Pi_00}) * ({epsilon_F}) / ({n})")
print("\nThe final numerical value is:")
print(result)
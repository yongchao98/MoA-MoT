import numpy as np

# Step 1: Define constants and convert units
M = 938.93  # MeV
Cs2 = 267.1
nu = 4
kF_fm = 1.42  # fm^-1
R_fm = 5.5    # fm
hbar_c = 197.3  # MeV fm

# Convert to MeV-based units
kF = kF_fm * hbar_c  # MeV
R = R_fm / hbar_c    # MeV^-1

# Step 2: Calculate M*
M_star = 0.78 * M

# Step 3: Calculate the momentum integral I_k
kF2 = kF**2
M_star2 = M_star**2
integrand_num = kF2 + 2 * M_star2
integrand_den = np.sqrt(kF2 + M_star2)
Ik = integrand_num / integrand_den - 2 * M_star

# Step 4: Estimate the NC parameter eta
# Assume eta * theta = 1 and sqrt(theta) = R (in fm)
# So, eta = 1 / R_fm^2 (in fm^-2)
# Convert eta to MeV^2
eta = (1 / R_fm**2) * (hbar_c**2)

# Step 5: Calculate the fractional change
# Fractional change = eta * (Cs2 / M**2) * nu * (3 * R * Ik) / (32 * pi**2)
coeff_eta = (Cs2 / M**2) * nu * (3 * R * Ik) / (32 * np.pi**2)
fractional_change = eta * coeff_eta

# Step 6: Convert to percentage change
percentage_change = fractional_change * 100

# Step 7: Print the final equation with numerical values
print("Calculation of the percentage change in effective nucleon mass due to NCQM.")
print("="*70)
print(f"Given constants:")
print(f"  M = {M:.2f} MeV")
print(f"  C_s^2 = {Cs2:.1f}")
print(f"  nu = {nu}")
print(f"  k_F = {kF_fm:.2f} fm^-1 = {kF:.2f} MeV")
print(f"  R = {R_fm:.1f} fm = {R:.5f} MeV^-1")
print("\nIntermediate calculations:")
print(f"  Standard effective mass, M* = 0.78 * M = {M_star:.2f} MeV")
print(f"  Momentum integral, Ik = [({kF:.2f}^2 + 2*{M_star:.2f}^2)/sqrt({kF:.2f}^2 + {M_star:.2f}^2)] - 2*{M_star:.2f} = {Ik:.4f} MeV")
print(f"  NC parameter assumption: eta*theta = 1, sqrt(theta) = R")
print(f"  Estimated eta = (hbar_c/R_fm)^2 = ({hbar_c:.1f}/{R_fm:.1f})^2 = {eta:.2f} MeV^2")
print("\nFinal Equation for Percentage Change:")
print("  % Change = (C_s^2 / M^2) * nu * (3 * R * Ik) / (32 * pi^2) * eta * 100")
print(f"  % Change = ({Cs2:.1f} / {M:.2f}^2) * {nu} * (3 * {R:.5f} * {Ik:.4f}) / (32 * pi^2) * {eta:.2f} * 100")

# Final calculated result
print("\nResult:")
print(f"The percentage change in the effective mass is: +{percentage_change:.4f}%")

# Choosing the closest answer from the options
options = {'A': 0.12, 'B': 5.0, 'C': 0.0005, 'D': -5.0, 'E': -0.12}
closest_option = min(options.items(), key=lambda item: abs(item[1] - percentage_change))

# The result is approximately +0.144%, which is closest to +0.12%.
print(f"This is closest to option {closest_option[0]} ({closest_option[1]:+}%).")
print("<<<A>>>")

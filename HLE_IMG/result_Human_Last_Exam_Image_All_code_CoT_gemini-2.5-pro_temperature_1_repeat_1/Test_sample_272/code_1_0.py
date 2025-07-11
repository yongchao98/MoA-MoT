import math

# Step 1: Define constants and parameters from the problem
# Physical mass of the Higgs boson in GeV
m_h_phys = 125.5
# Top quark Yukawa coupling
y_t = 0.95
# Current energy scale (Lambda_1) in TeV
lambda1_tev = 8.0
# Proposed new energy scale (Lambda_2) in PeV
lambda2_pev = 1.1

# Step 2: Convert all energy scales to a common unit (GeV) for consistency
# 1 TeV = 1,000 GeV
lambda1_gev = lambda1_tev * 1000
# 1 PeV = 1,000,000 GeV
lambda2_gev = lambda2_pev * 1000 * 1000

# Step 3: Formulate the calculation for the multiplicative factor
# The factor is the ratio of the bare Higgs mass at the two different scales.
# Factor = m_h_bare(Lambda_2) / m_h_bare(Lambda_1)
# Using m_h_bare^2 = m_h_phys^2 + (y_t^2 / (16*pi^2)) * Lambda^2

# Step 4: Calculate the components of the formula and display them
m_h_phys_sq = m_h_phys**2
k = y_t**2 / (16 * math.pi**2)

# Numerator term for the ratio's square (at Lambda_2)
numerator_term = m_h_phys_sq + k * lambda2_gev**2
# Denominator term for the ratio's square (at Lambda_1)
denominator_term = m_h_phys_sq + k * lambda1_gev**2

print("The formula for the multiplicative factor is:")
print("Factor = sqrt( (m_h_phys^2 + (y_t^2 / (16*pi^2)) * Lambda_2^2) / (m_h_phys^2 + (y_t^2 / (16*pi^2)) * Lambda_1^2) )")
print("\n--- Plugging in the numbers ---")
print(f"m_h_phys = {m_h_phys} GeV")
print(f"y_t = {y_t}")
print(f"Lambda_1 = {lambda1_gev} GeV")
print(f"Lambda_2 = {lambda2_gev} GeV")
print("-" * 30)

print("--- Calculating the components of the equation ---")
print(f"m_h_phys^2 = {m_h_phys_sq:.2f} GeV^2")
print(f"y_t^2 / (16*pi^2) = {k:.6f}")
print(f"Denominator = {m_h_phys_sq:.2f} + {k:.6f} * ({lambda1_gev})^2 = {denominator_term:.2f} GeV^2")
print(f"Numerator = {m_h_phys_sq:.2f} + {k:.6f} * ({lambda2_gev})^2 = {numerator_term:.2f} GeV^2")
print("-" * 30)


# Step 5: Calculate the final factor and round it to two decimal places
multiplicative_factor = math.sqrt(numerator_term / denominator_term)
rounded_factor = round(multiplicative_factor, 2)

print("--- Final Calculation ---")
print(f"Multiplicative factor = sqrt({numerator_term:.2f} / {denominator_term:.2f})")
print(f"The fine-tuning changes by a multiplicative factor of: {rounded_factor}")
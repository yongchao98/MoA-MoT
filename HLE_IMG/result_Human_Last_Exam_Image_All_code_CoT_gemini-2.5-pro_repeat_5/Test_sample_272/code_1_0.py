import math

# Plan: Calculate the multiplicative factor by which the fine-tuning measure changes.
# The factor is the ratio of the bare Higgs masses at the two energy scales.

# --- Given constants ---
m_h_phys = 125.5  # Physical Higgs mass in GeV
y_t = 0.95        # Top quark Yukawa coupling
Lambda_1 = 8 * 1e3  # Current energy scale (8 TeV) in GeV
Lambda_2 = 1.1 * 1e6 # Proposed energy scale (1.1 PeV) in GeV

# --- Calculation ---

# The relationship for the bare mass is:
# m_h_bare = sqrt(m_h_phys^2 + (y_t^2 / (16*pi^2)) * Lambda^2)

# The multiplicative factor is m_h_bare(Lambda_2) / m_h_bare(Lambda_1)

# Calculate the coefficient from the loop calculation
coeff = y_t**2 / (16 * math.pi**2)

# --- Calculate bare mass at Lambda_1 ---
m_h_bare_sq_1 = m_h_phys**2 + coeff * Lambda_1**2
m_h_bare_1 = math.sqrt(m_h_bare_sq_1)

# --- Calculate bare mass at Lambda_2 ---
m_h_bare_sq_2 = m_h_phys**2 + coeff * Lambda_2**2
m_h_bare_2 = math.sqrt(m_h_bare_sq_2)

# --- Calculate the multiplicative factor ---
factor = m_h_bare_2 / m_h_bare_1

# --- Output the results with intermediate steps ---
print("The multiplicative factor is the ratio of the bare Higgs mass at the new scale (Lambda_2) to the old scale (Lambda_1).")
print("\nStep 1: Calculate the bare Higgs mass at Lambda_1 = 8 TeV.")
print(f"m_h_bare(Lambda_1) = sqrt({m_h_phys:.1f}^2 + ({y_t:.2f}^2 / (16*pi^2)) * {Lambda_1:.0f}^2)")
print(f"m_h_bare(Lambda_1) = {m_h_bare_1:.2f} GeV")

print("\nStep 2: Calculate the bare Higgs mass at Lambda_2 = 1.1 PeV.")
print(f"m_h_bare(Lambda_2) = sqrt({m_h_phys:.1f}^2 + ({y_t:.2f}^2 / (16*pi^2)) * {Lambda_2:.0f}^2)")
print(f"m_h_bare(Lambda_2) = {m_h_bare_2:.2f} GeV")

print("\nStep 3: Calculate the multiplicative factor.")
print(f"Factor = m_h_bare(Lambda_2) / m_h_bare(Lambda_1)")
print(f"Factor = {m_h_bare_2:.2f} / {m_h_bare_1:.2f}")
print(f"The final multiplicative factor is: {factor:.2f}")
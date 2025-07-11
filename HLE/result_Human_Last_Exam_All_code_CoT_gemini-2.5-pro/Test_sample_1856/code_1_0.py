import math

# --- Given values ---
# Velocity of particle A in the lab frame
beta_A = 0.95

# --- Step 1: Calculate the Lorentz factor gamma ---
gamma_A = 1 / math.sqrt(1 - beta_A**2)

# --- Step 2: Define the momentum components of C in A's rest frame ---
# The four-momentum of C in the rest frame is (E_C*, P_C*/sqrt(2), 0, P_C*/sqrt(2)).
# We are given m_C << E_C, which implies C is ultra-relativistic, so E_C* is approximately P_C*.
# The momentum components are thus expressed as coefficients of the unknown P_C*.
# p_Cx* = (1/sqrt(2)) * P_C*
# p_z*  = (1/sqrt(2)) * P_C*
# E_C*  = 1 * P_C* (the approximation)

p_cx_star_coeff = 1 / math.sqrt(2)
p_z_star_coeff = 1 / math.sqrt(2)
E_c_star_coeff = 1.0 # Due to the E_C* ~ P_C* approximation

# --- Step 3: Apply Lorentz transformation to find momentum components in the lab frame ---
# Transverse momentum is invariant: p_Cx = p_Cx*
p_cx_lab_coeff = p_cx_star_coeff

# Longitudinal momentum transforms: p_Cz = gamma_A * (p_z* + beta_A * E_C*)
p_cz_lab_coeff = gamma_A * (p_z_star_coeff + beta_A * E_c_star_coeff)

# --- Step 4: Calculate the angle in the lab frame ---
# The angle is given by tan(theta) = p_Cx / p_Cz.
# Since P_C* is a common factor, it cancels out.
tan_theta = p_cx_lab_coeff / p_cz_lab_coeff
# Calculate the angle in radians
theta_rad = math.atan(tan_theta)
# Convert the angle to degrees
theta_deg = math.degrees(theta_rad)

# --- Step 5: Print the results and the equation ---
print("This script calculates the angle of a decay product in the lab frame.")
print(f"Given parent particle velocity beta_A = {beta_A}")
print(f"Calculated Lorentz factor gamma_A = {gamma_A:.4f}")
print("\nUsing the ultra-relativistic approximation E_C* = P_C*:")
print(f"The transverse momentum is p_Cx = {p_cx_lab_coeff:.4f} * P_C*")
print(f"The longitudinal momentum is p_Cz = {p_cz_lab_coeff:.4f} * P_C*")
print("\nThe final angle is calculated as:")
print(f"theta = arctan(p_Cx / p_Cz)")
print(f"theta = arctan({p_cx_lab_coeff:.4f} / {p_cz_lab_coeff:.4f})")
print(f"theta = {theta_deg:.3f} degrees")

# The final answer in the required format
# print(f"\n<<<{theta_deg:.3f}>>>")
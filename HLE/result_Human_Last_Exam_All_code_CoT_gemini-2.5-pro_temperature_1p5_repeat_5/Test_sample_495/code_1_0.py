import math

def calculate_energy_loss():
    """
    Calculates the energy loss per centimeter for α-particles in air.
    
    This calculation is based on Geiger's law for the range-energy relationship
    of α-particles (R = a * E^(3/2)).
    """
    
    # Initial given parameters
    E_initial = 8.5  # Initial energy in MeV
    R_total = 8.3    # Total range in cm
    x = 4.0          # Distance from the source in cm
    
    print("--- Given Parameters ---")
    print(f"Initial Energy (E_initial): {E_initial} MeV")
    print(f"Total Range (R_total): {R_total} cm")
    print(f"Distance from source (x): {x} cm")
    print("-" * 26)
    
    # Step 1: Calculate the remaining range of the α-particle after traveling x cm.
    R_remaining = R_total - x
    print("\n--- Step 1: Calculate Remaining Range ---")
    print(f"The remaining range at {x} cm is:")
    print(f"R_remaining = R_total - x = {R_total} - {x} = {R_remaining:.3f} cm")

    # Step 2: Calculate the energy of the α-particle at the distance x.
    # From R = a*E^(3/2), we can derive E = (R/a)^(2/3).
    # This leads to the ratio: E_at_x / E_initial = (R_remaining / R_total)^(2/3).
    E_at_x = E_initial * math.pow(R_remaining / R_total, 2.0/3.0)
    print("\n--- Step 2: Calculate Energy at 4 cm ---")
    print(f"The energy at {x} cm (E_at_x) is calculated as:")
    print(f"E_at_x = E_initial * (R_remaining / R_total)^(2/3)")
    print(f"E_at_x = {E_initial} * ({R_remaining:.3f} / {R_total})^(2/3) = {E_at_x:.3f} MeV")
    
    # Step 3: Calculate the energy loss per centimeter (|dE/dx|).
    # From Geiger's law, we can derive |dE/dx| = E / (1.5 * R).
    # At distance x, we use the energy E_at_x and the remaining range R_remaining.
    energy_loss_per_cm = E_at_x / (1.5 * R_remaining)
    print("\n--- Step 3: Calculate Energy Loss per cm ---")
    print("The energy loss per cm (|dE/dx|) is calculated using the energy and remaining range at that point:")
    print(f"|dE/dx| = E_at_x / (1.5 * R_remaining)")
    print(f"|dE/dx| = {E_at_x:.3f} MeV / (1.5 * {R_remaining:.3f} cm)")

    print("\n--- Final Result ---")
    print(f"The energy loss per centimetre for these α-particles at a distance of {x} cm from the source is: {energy_loss_per_cm:.3f} MeV/cm.")

calculate_energy_loss()
<<<0.853>>>
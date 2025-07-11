import math

def calculate_energy_loss():
    """
    Calculates the energy loss per centimeter for α-particles in standard air
    based on their initial energy and range, using Geiger's Law.
    """
    # --- Given values ---
    E0 = 8.5  # Initial energy in MeV
    R0 = 8.3  # Total range in cm
    x = 4.0   # Distance from the source in cm

    # --- Step 1: Calculate the energy at distance x ---
    # The remaining range (Rx) after traveling a distance x.
    Rx = R0 - x

    # The energy (Ex) at distance x, derived from Geiger's Law (R ∝ E^(3/2)).
    # (Rx / R0) = (Ex / E0)^(3/2)  =>  Ex = E0 * (Rx / R0)^(2/3)
    Ex = E0 * (Rx / R0)**(2/3)

    # --- Step 2: Calculate the energy loss rate |dE/dx| at distance x ---
    # By differentiating Geiger's Law, we get the rate of energy loss
    # as a function of the particle's energy E at that point.
    # |dE/dx| = (2/3) * (E0^(3/2) / R0) * E^(-1/2)
    # We evaluate this at E = Ex.
    energy_loss_rate = (2/3) * (E0**1.5 / R0) * (Ex**-0.5)

    # --- Print the results step-by-step ---
    print(f"Given initial energy (E₀) = {E0} MeV")
    print(f"Given total range (R₀) = {R0} cm")
    print(f"Calculating energy loss at distance (x) = {x} cm")
    print("-" * 40)

    print("Step 1: Find the particle's energy at 4 cm.")
    print(f"The remaining range (Rx) is R₀ - x = {R0} - {x} = {Rx:.3f} cm.")
    print(f"The energy at 4 cm (Ex) is E₀ * (Rx/R₀)^(2/3):")
    print(f"Ex = {E0:.3f} * ({Rx:.3f}/{R0:.3f})^(2/3) = {Ex:.3f} MeV")
    print("-" * 40)

    print("Step 2: Calculate the energy loss per cm (|dE/dx|) at 4 cm.")
    print("The formula for energy loss rate is |dE/dx| = (2/3) * (E₀^(1.5) / R₀) * Ex^(-0.5)")
    print("Plugging in the numbers:")
    print(f"|dE/dx| = (2/3) * ({E0:.3f}^(1.5) / {R0:.3f}) * {Ex:.3f}^(-0.5)")
    print(f"The calculated energy loss per centimetre is: {energy_loss_rate:.3f} MeV/cm")

calculate_energy_loss()
import math

def calculate_energy_loss():
    """
    Calculates the energy loss of α-particles in standard air based on Geiger's rule.
    """
    # Initial parameters
    E0 = 8.5      # Initial energy in MeV
    R_total = 8.3 # Total range in cm
    x = 4.0       # Distance from the source in cm

    # The formula for energy loss (-dE/dx) derived from Geiger's rule is:
    # -dE/dx = (2 * E₀) / (3 * R_total^(2/3) * (R_total - x)^(1/3))

    # Calculate the numerator of the formula
    numerator = 2 * E0

    # Calculate the remaining range
    R_rem = R_total - x
    
    # Calculate the denominator of the formula
    denominator = 3 * (R_total**(2/3)) * (R_rem**(1/3))

    # Calculate the final energy loss
    energy_loss = numerator / denominator

    # --- Output ---
    print("This script calculates the energy loss per centimetre for α-particles in air.")
    print("The calculation is based on the following formula derived from Geiger's rule:")
    print("-dE/dx = (2 * E₀) / (3 * R_total^(2/3) * (R_total - x)^(1/3))\n")

    print("Substituting the given values into the equation:")
    print(f"E₀ (Initial Energy) = {E0} MeV")
    print(f"R_total (Total Range) = {R_total} cm")
    print(f"x (Distance from source) = {x} cm\n")

    print("The equation becomes:")
    # This print statement shows the equation with the numbers plugged in, as requested.
    print(f"-dE/dx = (2 * {E0}) / (3 * {R_total}^(2/3) * ({R_total} - {x})^(1/3))")
    
    print(f"-dE/dx = {numerator} / {denominator:.4f}\n")

    print(f"The energy loss per centimetre for these α-particles at a distance of {x} cm is {energy_loss:.3f} MeV/cm.")

if __name__ == "__main__":
    calculate_energy_loss()
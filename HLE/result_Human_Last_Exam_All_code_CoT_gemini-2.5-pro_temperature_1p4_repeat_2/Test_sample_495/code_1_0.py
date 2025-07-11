import math

def calculate_energy_loss():
    """
    Calculates the energy loss per centimetre for α-particles in air
    based on the Geiger-Nuttall range-energy relationship.
    """
    # Given parameters
    E0 = 8.5  # Initial energy in MeV
    R0 = 8.3  # Range in cm
    x = 4.0   # Distance from the source in cm

    # Explanation of the method and formula
    print("The energy loss per centimetre (dE/dx) is calculated using a formula derived from the empirical Geiger's rule (R = k*E^(3/2)).")
    print("The derived formula for the magnitude of energy loss at a distance x is:")
    print("|dE/dx| = (2/3) * (E₀/R₀) * (1 - x/R₀)^(-1/3)")
    print("\n--------------------------------------------------")
    
    # Printing the calculation with substituted values
    print("Substituting the given values:")
    print(f"Initial Energy (E₀) = {E0} MeV")
    print(f"Total Range (R₀) = {R0} cm")
    print(f"Distance (x) = {x} cm")
    
    print("\nEquation:")
    # Using format specifiers to show the numbers in the equation clearly
    print(f"|dE/dx| = (2/3) * ({E0:.1f} / {R0:.1f}) / (1 - {x:.1f} / {R0:.1f})^(1/3)")
    
    # Perform the calculation
    try:
        base_term = 1 - (x / R0)
        if base_term <= 0:
            print("\nError: The distance x must be less than the total range R0.")
            return

        energy_loss_per_cm = (2/3) * (E0 / R0) / (base_term**(1/3))

        # Output the final result
        print("\n--------------------------------------------------")
        print("Result:")
        print(f"The energy loss per centimetre for these α-particles at a distance of {x:.1f} cm is: {energy_loss_per_cm:.3f} MeV/cm")

    except Exception as e:
        print(f"\nAn error occurred during calculation: {e}")

# Execute the function
calculate_energy_loss()
<<<0.850>>>
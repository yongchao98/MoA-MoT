import math

def simulate_melting_curves():
    """
    Simulates and prints melting curve data for two distinct populations
    and their bulk average to demonstrate the loss of heterogeneity information.
    """
    # --- Constants and Parameters ---
    R = 8.314  # Gas constant in J/(mol*K)

    # Population 1: A less stable duplex
    tm1_celsius = 65.0
    dh1 = 350000  # Enthalpy in J/mol

    # Population 2: A more stable duplex (e.g., different sequence or buffer)
    tm2_celsius = 75.0
    dh2 = 380000  # Enthalpy in J/mol

    print("--- Simulation Parameters ---")
    print(f"Population 1: Tm = {tm1_celsius}°C, dH = {dh1} J/mol")
    print(f"Population 2: Tm = {tm2_celsius}°C, dH = {dh2} J/mol")
    print("Bulk measurement is simulated as a 50/50 average of the two populations.\n")

    # Convert melting temperatures to Kelvin for the equation
    tm1_kelvin = tm1_celsius + 273.15
    tm2_kelvin = tm2_celsius + 273.15

    def calculate_fraction_unfolded(T_celsius, Tm_kelvin, dH):
        """
        Calculates the fraction of unfolded molecules at a given temperature
        using the van't Hoff equation for a two-state transition.
        """
        T_kelvin = T_celsius + 273.15
        # The equation for the equilibrium constant K is exp((dH/R) * (1/Tm - 1/T))
        # The fraction unfolded is 1 / (1 + K)
        exponent = (dH / R) * ((1 / Tm_kelvin) - (1 / T_kelvin))
        # Prevent math overflow for very large exponents
        if exponent > 700:
            return 0.0
        K = math.exp(exponent)
        return 1 / (1 + K)

    print("--- Simulated Melting Data ---")
    print(f"{'Temp (°C)':<12} | {'Pop 1 Unfolded':<18} | {'Pop 2 Unfolded':<18} | {'Bulk Unfolded':<15}")
    print("-" * 70)

    # Simulate the experiment over a range of temperatures
    for temp_c in range(55, 86, 2):
        frac_unfolded1 = calculate_fraction_unfolded(temp_c, tm1_kelvin, dh1)
        frac_unfolded2 = calculate_fraction_unfolded(temp_c, tm2_kelvin, dh2)
        
        # The "bulk" measurement is the average of the two populations
        bulk_unfolded = (frac_unfolded1 + frac_unfolded2) / 2.0

        print(f"{temp_c:<12.1f} | {frac_unfolded1:<18.4f} | {frac_unfolded2:<18.4f} | {bulk_unfolded:<15.4f}")

simulate_melting_curves()
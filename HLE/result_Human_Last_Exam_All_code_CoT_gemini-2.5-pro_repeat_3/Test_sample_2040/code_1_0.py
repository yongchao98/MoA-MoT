def calculate_degenerate_spdc(pump_wavelength_nm):
    """
    Calculates the wavelength of down-converted photons in degenerate SPDC.

    In degenerate Spontaneous Parametric Down-Conversion (SPDC), a pump photon
    splits into two identical photons (signal and idler). Energy conservation dictates
    that the energy of the pump photon is split equally between the two new photons.
    Since energy is inversely proportional to wavelength (E = hc/λ), the wavelength
    of the down-converted photons is double that of the pump photon.

    Equation: λ_signal = 2 * λ_pump

    Args:
        pump_wavelength_nm (float): The wavelength of the pump laser in nanometers.
    """
    # Calculate the wavelength of the degenerate signal/idler photons
    signal_wavelength_nm = 2 * pump_wavelength_nm
    idler_wavelength_nm = signal_wavelength_nm

    # Print the explanation and the final equation with values
    print("Spontaneous Parametric Down-Conversion (SPDC) Energy Conservation:")
    print("Principle: A high-energy 'pump' photon splits into two lower-energy 'signal' and 'idler' photons.")
    print("\nFor the degenerate case (signal and idler photons are identical):")
    print("The wavelength of the new photons is double the pump wavelength.")
    print("\nCalculation:")
    print(f"Given a pump wavelength (λ_pump) of {pump_wavelength_nm} nm...")

    # Output each number in the final equation
    print("\nThe final equation is: λ_signal = 2 * λ_pump")
    print(f"Result: {signal_wavelength_nm} nm = 2 * {pump_wavelength_nm} nm")

# --- User Input ---
# Define the pump laser wavelength in nanometers.
# A common value for SPDC experiments is 405 nm (blue-violet laser).
pump_wl = 405.0

# Run the calculation and print the results
calculate_degenerate_spdc(pump_wl)
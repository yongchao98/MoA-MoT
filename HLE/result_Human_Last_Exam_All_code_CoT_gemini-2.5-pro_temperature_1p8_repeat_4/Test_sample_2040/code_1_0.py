def calculate_spdc_wavelength():
    """
    Calculates the idler photon wavelength in SPDC based on energy conservation.

    The energy conservation law for SPDC is:
    1 / pump_wavelength = 1 / signal_wavelength + 1 / idler_wavelength

    This script rearranges the formula to solve for the idler wavelength:
    1 / idler_wavelength = 1 / pump_wavelength - 1 / signal_wavelength
    """

    # --- Inputs ---
    # Wavelengths are in nanometers (nm)
    pump_wavelength = 405
    signal_wavelength = 810

    # --- Calculation ---
    # Calculate the reciprocal of the idler wavelength
    inv_idler_wavelength = (1 / pump_wavelength) - (1 / signal_wavelength)

    # Calculate the final idler wavelength by taking the inverse
    idler_wavelength = 1 / inv_idler_wavelength

    # --- Output ---
    # Print the equation with all the numbers, as requested
    print("The SPDC energy conservation equation is: 1/λ_pump = 1/λ_signal + 1/λ_idler")
    print("Substituting the known values (in nm):")
    # Using the calculated idler wavelength to show the complete, solved equation
    print(f"1 / {pump_wavelength:.0f} = 1 / {signal_wavelength:.0f} + 1 / {idler_wavelength:.0f}")

    # Print the final answer
    print("\nBased on the calculation:")
    print(f"The wavelength of the idler photon is {idler_wavelength:.0f} nm.")


if __name__ == "__main__":
    calculate_spdc_wavelength()
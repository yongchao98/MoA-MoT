def calculate_spdc_degenerate_case(pump_wavelength_nm):
    """
    Demonstrates energy conservation in degenerate Spontaneous Parametric
    Downconversion (SPDC), a process not expected in centrosymmetric materials
    like boron nanosheets but fundamental to nonlinear optics.

    In SPDC, a high-energy pump photon splits into two lower-energy photons
    (called signal and idler). The conservation of energy is described by the equation:
    1 / pump_wavelength = 1 / signal_wavelength + 1 / idler_wavelength

    This function calculates the outcome for the 'degenerate' case, where the
    signal and idler photons have the same wavelength.

    Args:
        pump_wavelength_nm (float): The wavelength of the pump laser in nanometers.
    """

    # In the degenerate case, signal_wavelength equals idler_wavelength.
    # The equation simplifies to: 1 / pump = 1 / degenerate + 1 / degenerate
    # Which is: 1 / pump = 2 / degenerate
    # Therefore: degenerate_wavelength = 2 * pump_wavelength
    degenerate_wavelength_nm = 2 * pump_wavelength_nm

    # Display the results
    print(f"This script demonstrates energy conservation for a hypothetical SPDC process.")
    print(f"Given a pump photon with a wavelength of {pump_wavelength_nm} nm:")
    print(f"In degenerate SPDC, it would split into two photons, each with a wavelength of {degenerate_wavelength_nm} nm.")
    print("\n---")

    # Output the final equation with each number, as requested.
    # We use integer casting for a clean representation of the equation.
    pump = int(pump_wavelength_nm)
    degen = int(degenerate_wavelength_nm)
    
    print("The final energy conservation equation is:")
    print(f"1 / {pump} = 1 / {degen} + 1 / {degen}")
    print("---")

# We will use the wavelength of a common violet laser pointer as an example pump source.
PUMP_WAVELENGTH = 405.0

calculate_spdc_degenerate_case(PUMP_WAVELENGTH)
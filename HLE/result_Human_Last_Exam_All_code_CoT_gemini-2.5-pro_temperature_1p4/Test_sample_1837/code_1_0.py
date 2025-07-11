import numpy as np

def solve_pandora_mystery():
    """
    Verifies astrophysical measurements for the star Pandora using Planck's Law.

    The script checks the consistency of spectral radiance (B), wavelength (L),
    and temperature (T). If inconsistent, it identifies the most likely
    erroneous measurement and calculates its correct value.
    """
    # Physical constants
    h = 6.62607015e-34  # Planck constant (J*s)
    c = 299792458       # Speed of light (m/s)
    k_B = 1.380649e-23   # Boltzmann constant (J/K)

    # Given measurements
    l_nm = 400.0                # Wavelength in nanometers
    l_m = l_nm * 1e-9           # Wavelength in meters
    T_given = 9000.0            # Estimated surface temperature in Kelvin
    B_measured = 1.2e15         # Measured spectral radiance in W/(m^2 * sr * m)

    # 1. Calculate the theoretical spectral radiance B_calc for the given L and T
    try:
        exponent = (h * c) / (l_m * k_B * T_given)
        B_calc = (2 * h * c**2) / (l_m**5 * (np.exp(exponent) - 1))
    except OverflowError:
        B_calc = 0.0

    # 2. Compare B_calc with B_measured. Use a 5% tolerance for "looking ok".
    if np.isclose(B_calc, B_measured, rtol=0.05):
        print(0)
        return

    # 3. If inconsistent, determine the most likely error.
    # As reasoned in the plan, the estimated temperature is the most likely candidate.
    # We will now calculate the expected temperature T_calc if B and L were correct.
    
    # Rearrange Planck's law to solve for T
    # T = (h*c) / (lambda * k_B * ln( (2*h*c^2)/(B * lambda^5) + 1 ))
    try:
        log_argument = (2 * h * c**2) / (B_measured * l_m**5) + 1
        T_calc = (h * c) / (l_m * k_B * np.log(log_argument))
    except (ValueError, ZeroDivisionError):
        # This would happen if log_argument is not positive, which is physically unlikely
        # given the inputs, but good practice to handle.
        print("Could not calculate corrected temperature due to mathematical error.")
        return

    # 4. Output the conclusion and the verification equation
    print("The provided measurements are inconsistent with Planck's Law.")
    print("The surface temperature (T), which was an estimate, is the quantity most likely to have a measurement error.")
    
    T_final = int(round(T_calc))
    print(f"The expected correct value for T is {T_final} K.")
    
    print("\nVerification using Planck's Law with the corrected temperature:")
    print("B = (2 * h * c^2) / (lambda^5 * (exp((h * c) / (lambda * k_B * T)) - 1))")

    # Calculate the Right-Hand Side (RHS) of the equation with the new temperature
    # to show that it matches the measured spectral radiance B_measured.
    exponent_corrected = (h * c) / (l_m * k_B * T_calc)
    RHS = (2 * h * c**2) / (l_m**5 * (np.exp(exponent_corrected) - 1))

    # Print the full equation with all numbers plugged in
    print("Plugging in the values:")
    print(f"{B_measured:.3e} = (2 * {h:.3e} * ({c:.3e})^2) / (({l_m:.3e})^5 * (exp(({h:.3e} * {c:.3e}) / ({l_m:.3e} * {k_B:.3e} * {T_final})) - 1))")
    print(f"Evaluating the right side gives: {RHS:.3e} W/(m^2*sr*m), which now matches the measured value.")

    # Output the final answer in the format Xv
    print(f"\nT{T_final}")


solve_pandora_mystery()
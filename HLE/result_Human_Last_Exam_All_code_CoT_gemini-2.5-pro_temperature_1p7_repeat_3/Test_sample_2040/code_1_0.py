def solve_spdc_wavelengths():
    """
    Calculates the wavelengths of photons in spontaneous parametric downconversion (SPDC)
    based on the principle of energy conservation.
    """
    # Known values
    lambda_p = 405.0  # Wavelength of the pump laser in nm
    lambda_s_non_degenerate = 780.0 # Given signal wavelength for the non-degenerate case in nm

    print("Solving for photon wavelengths in Spontaneous Parametric Down-Conversion (SPDC).")
    print(f"Pump Wavelength (λ_pump): {lambda_p} nm")
    print("-" * 50)

    # Case 1: Degenerate SPDC
    # In this case, λ_signal = λ_idler = λ_degenerate.
    # The energy conservation equation 1/λ_pump = 1/λ_signal + 1/λ_idler becomes:
    # 1/λ_pump = 2/λ_degenerate  =>  λ_degenerate = 2 * λ_pump
    print("Case 1: Degenerate SPDC (photons have equal wavelengths)")
    lambda_degenerate = 2 * lambda_p
    
    # Output the equation with numbers
    print("The calculation is based on the equation: λ_degenerate = 2 * λ_pump")
    print(f"Final Equation: {lambda_degenerate:.1f} nm = 2 * {lambda_p:.1f} nm")
    print(f"Result: The wavelength of each generated photon is {lambda_degenerate:.1f} nm.")
    print("-" * 50)

    # Case 2: Non-degenerate SPDC
    # In this case, λ_signal and λ_idler are different.
    # We solve for λ_idler using: 1/λ_idler = 1/λ_pump - 1/λ_signal
    # λ_idler = 1 / (1/λ_pump - 1/λ_signal)
    print("Case 2: Non-degenerate SPDC")
    print(f"Given Signal Wavelength (λ_signal): {lambda_s_non_degenerate} nm")
    
    # Perform calculation
    inv_lambda_p = 1 / lambda_p
    inv_lambda_s = 1 / lambda_s_non_degenerate
    inv_lambda_i = inv_lambda_p - inv_lambda_s
    lambda_idler = 1 / inv_lambda_i

    # Output the equation with numbers
    print("The calculation is based on the equation: 1/λ_idler = 1/λ_pump - 1/λ_signal")
    print(f"Final Equation: 1/{lambda_idler:.1f} nm = 1/{lambda_p:.1f} nm - 1/{lambda_s_non_degenerate:.1f} nm")
    print(f"Result: The wavelength of the idler photon is {lambda_idler:.1f} nm.")
    print("-" * 50)

# Execute the function
solve_spdc_wavelengths()
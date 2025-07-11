import sys

def calculate_spdc_idler_wavelength(lambda_pump, lambda_signal):
    """
    Calculates the idler wavelength in Spontaneous Parametric Down-Conversion (SPDC)
    based on energy conservation.

    The relationship is: 1/lambda_pump = 1/lambda_signal + 1/lambda_idler

    Args:
        lambda_pump (float): Wavelength of the pump photon in nanometers (nm).
        lambda_signal (float): Wavelength of the signal photon in nanometers (nm).
    """

    # The energy of the signal photon cannot be greater than the pump photon,
    # which means its wavelength cannot be shorter.
    if lambda_signal < lambda_pump:
        print(f"Error: Signal wavelength ({lambda_signal} nm) cannot be shorter than pump wavelength ({lambda_pump} nm).", file=sys.stderr)
        return

    # Rearrange the formula to solve for the idler wavelength:
    # 1/lambda_idler = 1/lambda_pump - 1/lambda_signal
    try:
        inv_lambda_pump = 1 / lambda_pump
        inv_lambda_signal = 1 / lambda_signal
        inv_lambda_idler = inv_lambda_pump - inv_lambda_signal

        # If the result is zero or negative, it means the wavelengths are equal or invalid.
        if inv_lambda_idler <= 0:
             print(f"Error: Invalid calculation. Check if signal wavelength is equal to pump wavelength.", file=sys.stderr)
             return

        lambda_idler = 1 / inv_lambda_idler

        # Output the results
        print(f"Pump Wavelength: {lambda_pump} nm")
        print(f"Signal Wavelength: {lambda_signal} nm")
        print("-" * 30)
        print(f"Calculated Idler Wavelength: {lambda_idler:.2f} nm")
        print("-" * 30)
        print("Verifying the energy conservation equation:")
        # Output each number in the final equation
        print(f"1 / {lambda_pump} = 1 / {lambda_signal} + 1 / {lambda_idler:.2f}")
        print(f"{inv_lambda_pump:.6f} = {inv_lambda_signal:.6f} + {inv_lambda_idler:.6f}")

    except ZeroDivisionError:
        print("Error: Wavelength values cannot be zero.", file=sys.stderr)


if __name__ == '__main__':
    # --- Example Values ---
    # A common green laser pump photon
    pump_wavelength_nm = 532.0
    # A chosen near-infrared signal photon
    signal_wavelength_nm = 900.0

    calculate_spdc_idler_wavelength(pump_wavelength_nm, signal_wavelength_nm)

import numpy as np

# Step 1: Define the known wavelengths in nanometers (nm).
# Let's assume a green pump laser and a near-infrared signal photon.
lambda_pump = 532.0
lambda_signal = 810.0

print(f"Pump Wavelength (λ_pump): {lambda_pump} nm")
print(f"Signal Wavelength (λ_signal): {lambda_signal} nm")
print("-" * 30)

# Step 2: Calculate the idler wavelength using the energy conservation formula.
# The formula is: 1/λ_pump = 1/λ_signal + 1/λ_idler
# Rearranging for λ_idler: λ_idler = 1 / (1/λ_pump - 1/λ_signal)

# First, check if the signal wavelength is physically possible (must be longer than the pump wavelength).
if lambda_signal <= lambda_pump:
    print("Error: Signal wavelength must be longer than the pump wavelength for SPDC.")
else:
    # Calculate the terms for the inverse wavelengths.
    inv_lambda_pump = 1 / lambda_pump
    inv_lambda_signal = 1 / lambda_signal

    # Calculate the inverse idler wavelength.
    inv_lambda_idler = inv_lambda_pump - inv_lambda_signal

    # Calculate the final idler wavelength.
    lambda_idler = 1 / inv_lambda_idler

    print(f"Calculated Idler Wavelength (λ_idler): {lambda_idler:.2f} nm\n")

    # As requested, printing the final equation with all the numbers.
    print("Energy Conservation Equation Check:")
    print(f"1 / {lambda_pump:.2f} = 1 / {lambda_signal:.2f} + 1 / {lambda_idler:.2f}")
    print(f"{inv_lambda_pump:.6f} = {inv_lambda_signal:.6f} + {inv_lambda_idler:.6f}")
    print(f"{inv_lambda_pump:.6f} = {inv_lambda_signal + inv_lambda_idler:.6f}")

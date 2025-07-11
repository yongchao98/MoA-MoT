import sys

# Plan: Demonstrate the energy conservation principle in Spontaneous Parametric Down-Conversion (SPDC).
# 1. Define the wavelength of the incoming high-energy photon (the "pump").
# 2. Define the wavelength of one of the down-converted photons (the "signal").
# 3. Calculate the wavelength of the other down-converted photon (the "idler") based on the law of energy conservation.
#    The equation is: 1/lambda_pump = 1/lambda_signal + 1/lambda_idler
# 4. Print the full equation with the calculated numbers.

# Step 1 & 2: Define pump and signal wavelengths in nanometers (nm).
# A common green laser pointer has a wavelength of 532 nm.
pump_wavelength = 532.0
# Let's assume we detect a signal photon in the near-infrared at 810 nm.
signal_wavelength = 810.0

# The energy of a photon is inversely proportional to its wavelength.
# For SPDC to be possible, the pump photon must have more energy (i.e., a shorter wavelength)
# than the signal photon.
if pump_wavelength >= signal_wavelength:
    print(f"Error: The pump wavelength ({pump_wavelength} nm) must be shorter than the signal wavelength ({signal_wavelength} nm).")
    sys.exit()

# Step 3: Calculate the idler wavelength.
# Rearranging the formula: 1/lambda_idler = 1/lambda_pump - 1/lambda_signal
inv_pump = 1 / pump_wavelength
inv_signal = 1 / signal_wavelength
inv_idler = inv_pump - inv_signal

# The idler wavelength is the reciprocal of this result.
idler_wavelength = 1 / inv_idler

# Step 4: Print the final equation with all the numbers.
print("Demonstration of Energy Conservation in SPDC:")
print("Equation: 1/λ_pump = 1/λ_signal + 1/λ_idler")
print("\nGiven values:")
print(f"Pump Wavelength (λ_pump): {pump_wavelength} nm")
print(f"Signal Wavelength (λ_signal): {signal_wavelength} nm")
print("\nCalculation:")
print(f"1 / {pump_wavelength} = 1 / {signal_wavelength} + 1 / λ_idler")
print(f"1 / λ_idler = {inv_pump:.6f} - {inv_signal:.6f}")
print(f"1 / λ_idler = {inv_idler:.6f}")
print("\nFinal Result:")
# We format the numbers back into the original equation for the final output.
print(f"1 / {pump_wavelength} nm = 1 / {signal_wavelength} nm + 1 / {idler_wavelength:.2f} nm")

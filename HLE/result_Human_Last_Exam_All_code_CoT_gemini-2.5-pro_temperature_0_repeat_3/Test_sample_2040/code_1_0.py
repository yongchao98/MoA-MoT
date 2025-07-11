# Plan: Calculate the idler wavelength in a Spontaneous Parametric Down-Conversion (SPDC) process.
# The calculation is based on the law of energy conservation.

# 1. Define the known wavelengths in nanometers (nm).
# We'll use a common green laser pump and a near-infrared signal wavelength.
pump_wavelength = 532.0
signal_wavelength = 800.0

# 2. State the physical principle.
# The energy of a photon is inversely proportional to its wavelength (E = hc/λ).
# Energy conservation in SPDC means: E_pump = E_signal + E_idler
# This translates to the wavelength equation: 1/λ_pump = 1/λ_signal + 1/λ_idler.

# 3. Rearrange the equation to solve for the unknown idler wavelength.
# 1/λ_idler = 1/λ_pump - 1/λ_signal

# 4. Perform the calculation.
# Calculate the reciprocal of the idler wavelength first.
reciprocal_idler_wavelength = (1 / pump_wavelength) - (1 / signal_wavelength)

# Now, calculate the idler wavelength by taking the inverse.
idler_wavelength = 1 / reciprocal_idler_wavelength

# 5. Print the final equation with all the numbers, as requested.
print("The energy conservation equation for this SPDC process is:")
# The f-string below substitutes the variable values into the string.
# The idler wavelength is formatted to one decimal place for clarity.
print(f"1 / {pump_wavelength} = 1 / {signal_wavelength} + 1 / {idler_wavelength:.1f}")

# 6. Print the result separately for the user.
print(f"\nGiven a pump wavelength of {pump_wavelength} nm and a signal wavelength of {signal_wavelength} nm, the calculated idler wavelength is {idler_wavelength:.1f} nm.")

# This script demonstrates the principle of broadband CARS microscopy.
# It shows that different vibrational modes in a sample produce
# distinguishable frequencies in the anti-Stokes signal.

# Frequencies are in wavenumbers (cm^-1), a common unit in spectroscopy.

# Let's define the frequency of the probe beam. In many CARS setups,
# the probe and pump beams are from the same laser source.
omega_probe = 20000  # cm^-1

# Now, let's consider two different molecular vibrations present in a sample.
# These represent the 'information' we want to distinguish.
# The vibrational frequency is represented by Omega (Ω).
vibration_1_Omega = 2850  # A typical C-H chemical bond stretch vibration
vibration_2_Omega = 3000  # A different C-H or aromatic ring vibration

print("Demonstrating the Broadband CARS Principle:")
print("The anti-Stokes signal (ω_as) is generated based on the sample's vibration (Ω) and the probe beam (ω_probe).")
print("The relationship is: ω_as = ω_probe + Ω\n")

# --- Case 1: Detecting the First Vibration ---
# In broadband CARS, a range of vibrations are excited simultaneously.
# Let's see what signal is produced by the first vibration.
omega_as_1 = omega_probe + vibration_1_Omega

print("Case 1: Detecting the first vibration")
print(f"For a vibration Ω_1 = {vibration_1_Omega} cm^-1:")
# We print the final equation with the specific numbers
print(f"ω_as_1 = {omega_probe} + {vibration_1_Omega} = {omega_as_1} cm^-1\n")


# --- Case 2: Detecting the Second Vibration ---
# The broadband excitation also covers the second vibrational mode.
# Let's calculate its resulting anti-Stokes frequency.
omega_as_2 = omega_probe + vibration_2_Omega

print("Case 2: Detecting the second vibration")
print(f"For a vibration Ω_2 = {vibration_2_Omega} cm^-1:")
# We print the final equation with the specific numbers
print(f"ω_as_2 = {omega_probe} + {vibration_2_Omega} = {omega_as_2} cm^-1\n")


print("--- Conclusion ---")
print(f"The two different vibrations, {vibration_1_Omega} cm^-1 and {vibration_2_Omega} cm^-1, produced two different anti-Stokes signals at {omega_as_1} cm^-1 and {omega_as_2} cm^-1.")
print("When collected by a spectrometer, these signals are clearly separate.")
print("Therefore, the generated anti-Stokes beam contains distinguishable vibrational information.")
<<<C>>>
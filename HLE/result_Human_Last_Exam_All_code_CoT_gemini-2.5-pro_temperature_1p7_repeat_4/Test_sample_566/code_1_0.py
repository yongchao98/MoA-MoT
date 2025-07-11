# All frequencies are in wavenumbers (cm^-1) for convenience.

# Define the frequency of a narrowband Stokes beam.
# (e.g., from a 1064 nm Nd:YAG laser)
w_s = 9398 

# Define two different molecular vibrational frequencies (Raman shifts) in the sample.
Omega_1 = 2900  # A typical C-H bond stretch
Omega_2 = 1650  # A typical C=C bond stretch

# 1. Determine which frequency component from the BROADBAND PUMP excites each vibration.
# The resonance condition is Omega = w_p - w_s  =>  w_p = Omega + w_s
w_p1 = Omega_1 + w_s
w_p2 = Omega_2 + w_s

print("A broadband pump source provides the necessary frequencies to excite multiple vibrations.")
print(f"To excite vibration Omega_1 ({Omega_1} cm^-1), the required pump frequency is w_p1 = {w_p1} cm^-1.")
print(f"To excite vibration Omega_2 ({Omega_2} cm^-1), the required pump frequency is w_p2 = {w_p2} cm^-1.")
print("-" * 50)

# 2. Calculate the frequency of the resulting anti-Stokes signal for each vibration.
# The signal frequency is w_as = 2*w_p - w_s (assuming pump and probe are the same)
w_as1 = 2 * w_p1 - w_s
w_as2 = 2 * w_p2 - w_s

print("Each vibration generates an anti-Stokes signal at a unique, distinguishable frequency.")
print("\nFinal equation for the signal from vibration 1:")
print(f"w_as1 = 2 * {w_p1} - {w_s} = {w_as1} cm^-1")

print("\nFinal equation for the signal from vibration 2:")
print(f"w_as2 = 2 * {w_p2} - {w_s} = {w_as2} cm^-1")
print("-" * 50)

# 3. Verify that the information is distinguishable.
if w_as1 != w_as2:
    print(f"\nResult: The signal frequencies {w_as1} cm^-1 and {w_as2} cm^-1 are different.")
    print("Therefore, the broadband anti-Stokes beam contains distinguishable vibrational information.")
else:
    print("\nResult: The signal frequencies are the same, information is not distinguishable.")

import math

# --- Parameters ---
# Let's model a simplified coil for a 1.5T MRI scanner.
# Operating frequency should be around 63.87 MHz.

# L: Inductance of the coil loop (in Henries)
L = 100e-9  # 100 nH

# C: Tuning capacitance for the operating frequency (in Farads)
# We can calculate C to resonate with L at ~63.87 MHz
# f = 1 / (2*pi*sqrt(L*C)) => C = 1 / ((2*pi*f)^2 * L)
C = 62e-12  # 62 pF

# C_detune: An additional capacitor switched in by a PIN diode
# to detune the coil.
C_detune = 500e-12 # 500 pF

pi = math.pi

# --- Calculations ---

# 1. Tuned State (Receive Mode)
# In this mode, the detuning circuit is off, and the coil resonates at its
# intended operating frequency. The total capacitance is just C.
# The equation is: f_tuned = 1 / (2 * pi * sqrt(L * C))
f_tuned = 1 / (2 * pi * math.sqrt(L * C))


# 2. Detuned State (Transmit/Safe Mode)
# The PIN diode is activated, adding C_detune in parallel with C.
# The total capacitance is now C + C_detune.
# The equation is: f_detuned = 1 / (2 * pi * sqrt(L * (C + C_detune)))
C_total_detuned = C + C_detune
f_detuned = 1 / (2 * pi * math.sqrt(L * C_total_detuned))

# --- Output the Results ---
# The output clearly shows each number used in the final equations.

print("This model demonstrates how a coil's resonant frequency is shifted by an active detuning circuit.")
print("-" * 70)

# Print the tuned equation and result
print("Tuned State (Receive Mode):")
print(f"f_tuned = 1 / (2 * pi * sqrt(L * C))")
print(f"f_tuned = 1 / (2 * {pi:.4f} * sqrt({L:.3e} H * {C:.3e} F))")
print(f"Result -> The coil's operating frequency is {f_tuned / 1e6:.2f} MHz.\n")


# Print the detuned equation and result
print("Detuned State (Safe Mode, as measured on the bench without power):")
print(f"f_detuned = 1 / (2 * pi * sqrt(L * (C + C_detune)))")
print(f"f_detuned = 1 / (2 * {pi:.4f} * sqrt({L:.3e} H * ({C:.3e} F + {C_detune:.3e} F)))")
print(f"Result -> The detuned frequency is {f_detuned / 1e6:.2f} MHz.")

print("-" * 70)
print(f"\nConclusion: If you are looking for a resonance at {f_tuned / 1e6:.2f} MHz,")
print(f"you won't see it because the unpowered coil is actually resonant at {f_detuned / 1e6:.2f} MHz.")

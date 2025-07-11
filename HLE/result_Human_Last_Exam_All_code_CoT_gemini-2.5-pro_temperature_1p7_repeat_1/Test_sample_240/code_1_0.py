import numpy as np

def calculate_resonant_frequency(L, C):
    """Calculates the resonant frequency of an LC circuit."""
    return 1 / (2 * np.pi * np.sqrt(L * C))

def calculate_capacitance(f, L):
    """Calculates the required capacitance for a given resonant frequency and inductance."""
    return 1 / ((2 * np.pi * f)**2 * L)

# --- Step 1: Define parameters for a typical 1.5T MRI receive coil ---
# The Larmor frequency for protons at 1.5 Tesla is approximately 63.87 MHz.
f_operating = 63.87e6  # Operating frequency in Hz
# Let's assume a typical inductance for a small surface coil.
L_coil = 200e-9       # Coil inductance in Henries (200 nH)

print(f"--- Calculating for a Tuned Coil (Receive Mode) ---")
print(f"Target Operating Frequency (f): {f_operating / 1e6:.2f} MHz")
print(f"Assumed Coil Inductance (L): {L_coil * 1e9:.1f} nH")

# --- Step 2: Calculate the required capacitance for the coil to be 'tuned' ---
C_tuned = calculate_capacitance(f_operating, L_coil)
print(f"\nEquation: C = 1 / ((2 * pi * f)^2 * L)")
print(f"Calculation: C = 1 / ((2 * pi * {f_operating:.2e})^2 * {L_coil:.2e})")
print(f"Result: The required 'tuned' capacitance (C_tuned) is {C_tuned * 1e12:.2f} pF.")

# --- Step 3: Verify the resonant frequency with the calculated values ---
f_calculated = calculate_resonant_frequency(L_coil, C_tuned)
print(f"\nVerifying the resonant frequency with these values:")
print(f"Equation: f = 1 / (2 * pi * sqrt(L * C))")
print(f"Calculation: f = 1 / (2 * pi * sqrt({L_coil:.2e} * {C_tuned:.2e}))")
print(f"Result: The calculated resonant frequency is {f_calculated / 1e6:.2f} MHz, which matches the target.")

print("\n" + "="*50 + "\n")

# --- Step 4: Simulate a 'detuned' state ---
# A decoupling circuit might, for example, short out a main tuning capacitor or
# add another capacitor in parallel, significantly changing the total capacitance.
# Let's simulate adding 10 pF from a PIN diode circuit.
C_decoupling = 10e-12 # 10 pF
C_detuned = C_tuned + C_decoupling

print(f"--- Calculating for a Detuned Coil (Transmit or Inactive Mode) ---")
print(f"A decoupling circuit activates, changing the total capacitance.")
print(f"New Capacitance (C_detuned): {C_detuned * 1e12:.2f} pF (Original {C_tuned * 1e12:.2f} pF + Decoupling {C_decoupling * 1e12:.2f} pF)")

# --- Step 5: Calculate the new, detuned resonant frequency ---
f_detuned = calculate_resonant_frequency(L_coil, C_detuned)
print(f"\nCalculating the new, detuned resonant frequency:")
print(f"Equation: f = 1 / (2 * pi * sqrt(L * C))")
print(f"Calculation: f = 1 / (2 * pi * sqrt({L_coil:.2e} * {C_detuned:.2e}))")
print(f"Result: The new resonant frequency is {f_detuned / 1e6:.2f} MHz.")
print(f"\nThis new frequency is far from the operating frequency of {f_operating / 1e6:.2f} MHz.")
print("This is why the resonance peak seems to 'disappear' on your VNA.")

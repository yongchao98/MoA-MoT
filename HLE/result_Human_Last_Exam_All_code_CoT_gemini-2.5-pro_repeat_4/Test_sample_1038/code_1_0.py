import sys

# This script calculates the maximum allowable gate leakage current to meet the design specification.
# The specification requires that the voltage on the gate capacitor does not change
# by more than 1% per second due to leakage.

# --- Parameters ---
# Gate capacitance in Farads. We use 1pF (1e-12 F) as a representative value from option A.
C_gate = 1.0e-12

# Assumed effective voltage swing stored on the capacitor in Volts.
# A 0.5V swing is a reasonable assumption for a 1.2V system.
V_swing = 0.5

# Allowed leakage percentage per second.
leakage_rate_percent = 1.0

# --- Calculation ---
# 1. Calculate the allowed percentage as a decimal fraction.
leakage_rate_decimal = leakage_rate_percent / 100.0

# 2. Calculate the maximum allowed voltage drop per second (dV/dt).
# dV/dt = V_swing * leakage_rate
# The time 'dt' is 1 second by definition of the spec.
dV_per_second = V_swing * leakage_rate_decimal

# 3. Calculate the maximum leakage current using the capacitor current formula: I = C * dV/dt.
I_leak_max_amps = C_gate * dV_per_second

# 4. Convert the result to femtoamperes (fA) for easier interpretation, as 1 fA = 1e-15 A.
I_leak_max_femtoamps = I_leak_max_amps / 1.0e-15

# --- Output ---
# The final output presents the equation with the values plugged in.
# sys.stdout.write is used for cleaner output without extra newlines from print().
sys.stdout.write("Calculating the maximum leakage current (I_leak) allowed:\n")
sys.stdout.write("I_leak = C_gate * (V_swing * leakage_rate / 1s)\n")
sys.stdout.write(f"I_leak = {C_gate:.1e} F * ({V_swing:.1f} V * {leakage_rate_percent / 100:.2f} / 1s)\n")
sys.stdout.write(f"I_leak = {C_gate:.1e} F * {dV_per_second:.3f} V/s\n")
sys.stdout.write(f"I_leak = {I_leak_max_amps:.1e} A\n\n")
sys.stdout.write(f"This means the maximum allowed leakage current is {I_leak_max_femtoamps:.1f} femtoamperes (fA).\n")
sys.stdout.write("This extremely low value highlights the difficulty of meeting the leakage specification.\n")

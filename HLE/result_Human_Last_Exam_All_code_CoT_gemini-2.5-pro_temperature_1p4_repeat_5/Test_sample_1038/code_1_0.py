import math

# --- Design Constraints ---
V_SUPPLY = 1.2  # volts
VT = 0.45  # volts
RESET_TIME_TARGET = 5e-6  # 5 microseconds
LEAKAGE_CURRENT = 10e-15 # Assume 10 fA leakage current
LEAKAGE_TOLERANCE_PERCENT = 1.0 # % per second drift allowed
# Voltage on the cap relative to source to set subthreshold bias
V_BIAS = 0.3 # Assume a 300mV gate-source bias is needed

def evaluate_strategy(name, R_on, C_gate):
    """
    Calculates reset time and leakage drift for a given R_on and C_gate.
    Assuming reset is charging the cap to V_BIAS, we can approximate the
    time with a few time constants. Let's use 5*tau for settling.
    t_reset = 5 * R_on * C_gate
    """
    t_reset = 5 * R_on * C_gate
    
    # Voltage drift in one second
    v_drift_per_sec = LEAKAGE_CURRENT * 1.0 / C_gate
    
    # Drift as a percentage of the bias voltage
    percent_drift_per_sec = (v_drift_per_sec / V_BIAS) * 100
    
    # Check if constraints are met
    reset_ok = t_reset <= RESET_TIME_TARGET
    leakage_ok = percent_drift_per_sec <= LEAKAGE_TOLERANCE_PERCENT
    
    print(f"--- Strategy: {name} ---")
    print(f"Assumed On-Resistance (R_on): {R_on:,.0f} Ohms")
    print(f"Assumed Gate Capacitance (C_gate): {C_gate*1e12:.2f} pF")
    print("\nResults:")
    print(f"Calculated Reset Time: {t_reset*1e6:.2f} us (Target: < {RESET_TIME_TARGET*1e6:.2f} us) -> Met: {reset_ok}")
    print(f"Calculated Voltage Drift: {percent_drift_per_sec:.4f} % per second (Target: < {LEAKAGE_TOLERANCE_PERCENT:.1f} %) -> Met: {leakage_ok}")
    print(f"Overall Balance Achieved: {reset_ok and leakage_ok}")
    print("-" * 25 + "\n")

# Strategy A: Standard transistor (moderate R_on), small capacitor to meet reset time
# A standard reset transistor might have an R_on of 1 kOhm
R_on_A = 1000  # Ohms
C_gate_A = 1e-12   # 1 pF
evaluate_strategy("A (Small C_gate)", R_on_A, C_gate_A)

# Strategy E: Split-gate transistor (very low R_on), large capacitor for leakage tolerance
# The parallel structure could halve the R_on or better.
R_on_E = 500   # Ohms
# We can now use a larger capacitor because R_on is lower
# Let's find a C that just meets the reset time requirement with this new R_on
# 5e-6 = 5 * 500 * C_gate  => C_gate = 5e-6 / 2500 = 2e-9 F = 2000 pF
# Let's use a practical value like 10 pF, which is much larger than A's cap
C_gate_E = 10e-12  # 10 pF
evaluate_strategy("E (Split-Gate, Larger C_gate)", R_on_E, C_gate_E)

print("Analysis:")
print("Strategy A prioritizes a fast reset with a small capacitor (1.00 pF), but fails the leakage test dramatically.")
print("Strategy E uses a split-gate design to achieve a lower on-resistance during reset (e.g., 500 Ohms).")
print("This lower resistance allows for a much larger capacitor (10.00 pF) to be used while still easily meeting the reset time target.")
print("The larger capacitor in Strategy E makes the circuit far more robust against leakage, successfully balancing the conflicting requirements.")

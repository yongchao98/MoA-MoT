import numpy as np

# --- Problem Parameters (in mV and ms) ---
tau_m = 20.0      # Membrane time constant (ms)
tau_ref = 2.0     # Refractory period (ms)
V_th = 20.0       # Voltage threshold (mV)
V_reset = 10.0    # Voltage reset (mV)
J_E = 0.1         # Excitatory synaptic efficacy (mV)
K_E = 1000        # Number of excitatory connections
prop_I_to_E = 0.25# Proportion of inhibitory connections
g = 4.0           # Relative inhibition to excitation
mu_ext = 30.0     # External input mean depolarization (mV)

# --- Derived Parameters ---
K_I = K_E * prop_I_to_E  # Number of inhibitory connections (250)
J_I = g * J_E          # Inhibitory synaptic efficacy (0.4 mV)

# --- Step 1: Calculate the Mean Membrane Potential (mu) ---
# The mean depolarization from recurrent connections is proportional to:
# K_E * J_E - K_I * J_I
# = 1000 * 0.1 - 250 * 0.4
# = 100 - 100 = 0
# Because this term is zero, the recurrent input does not affect the mean potential (a "balanced network").
# Thus, the total mean input is just the external input.
mu = mu_ext

# --- Step 2: Calculate the Firing Rate ---
# Since mu (30 mV) > V_th (20 mV), the neuron is in the suprathreshold regime.
# We use the formula for the time it takes for the voltage to drift from V_reset to V_th.
# T_drift = tau_m * ln((mu - V_reset) / (mu - V_th))

# Calculate the terms inside the logarithm
mu_minus_V_reset = mu - V_reset
mu_minus_V_th = mu - V_th
ratio = mu_minus_V_reset / mu_minus_V_th

# Calculate drift time in ms
T_drift = tau_m * np.log(ratio)

# The total period is the drift time plus the refractory period (in ms)
total_period_ms = tau_ref + T_drift

# The firing rate in Hz is 1000 / total_period_ms
firing_rate_Hz = 1000.0 / total_period_ms

# --- Step 3: Print the Final Equation and Result ---
# The prompt requires printing each number in the final equation.

print("Final Equation for Firing Rate:")
print(f"rate (Hz) = 1000 / (τ_ref + τ_m * ln((μ - V_reset) / (μ - V_th)))")
print(f"rate (Hz) = 1000 / ({tau_ref:.1f} + {tau_m:.1f} * ln(({mu:.1f} - {V_reset:.1f}) / ({mu:.1f} - {V_th:.1f})))")
print(f"rate (Hz) = 1000 / ({tau_ref:.1f} + {tau_m:.1f} * ln({mu_minus_V_reset:.1f} / {mu_minus_V_th:.1f}))")
print(f"rate (Hz) = 1000 / ({tau_ref:.1f} + {tau_m:.1f} * ln({ratio:.1f}))")
print(f"rate (Hz) = 1000 / ({tau_ref:.1f} + {T_drift:.2f})")
print(f"rate (Hz) = 1000 / {total_period_ms:.2f}")
print(f"rate (Hz) = {firing_rate_Hz:.2f}")

# The problem asks for the answer as an integer.
final_answer = int(round(firing_rate_Hz))
print(f"\nThe final firing rate, rounded to the nearest integer, is: {final_answer} Hz")

<<<63>>>
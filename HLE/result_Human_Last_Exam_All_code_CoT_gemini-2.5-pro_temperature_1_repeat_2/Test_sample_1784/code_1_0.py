import math

# --- Step 1: Define Parameters in standard units (Volts, Seconds) ---
tau = 20e-3         # Membrane time constant (s)
V_reset = 10e-3     # Voltage reset (V)
V_th = 20e-3        # Voltage threshold (V)
tau_ref = 2e-3      # Refractory period (s)
V_ext = 30e-3       # External input resulting in this steady-state potential (V)
J_E = 0.1e-3        # Excitatory synaptic efficacy (V)
K_E = 1000          # Number of excitatory connections
K_I_prop = 0.25     # Proportion of inhibitory connections relative to K_E
g = 4.0             # Relative strength of inhibition to excitation (J_I / J_E)

# --- Step 2: Determine Mean Membrane Potential (μ) ---
# Calculate derived parameters
K_I = K_I_prop * K_E
J_I = g * J_E

# In a self-consistent state where excitatory (r_E) and inhibitory (r_I) neurons fire at the same rate r,
# the mean potential would be: μ = τ * r * (K_E * J_E - K_I * J_I) + V_ext
# Let's check the term in the parenthesis, which represents the net recurrent input.
net_recurrent_drive = K_E * J_E - K_I * J_I

# Because the net recurrent drive is 0 (1000 * 0.1e-3 - 250 * 0.4e-3 = 0.1e-3 - 0.1e-3 = 0),
# the network is in a balanced state. The average membrane potential is thus determined by the external input.
mu = V_ext

# --- Step 3: Apply Firing Rate Formula and Calculate ---
# The firing rate 'r' is given by: r = 1 / (τ_ref + τ * log((μ - V_reset) / (μ - V_th)))
# We must ensure μ > V_th for the neuron to fire.
# Here, μ = 30 mV and V_th = 20 mV, so the condition is met.

log_argument = (mu - V_reset) / (mu - V_th)
natural_log = math.log(log_argument)
denominator = tau_ref + tau * natural_log
firing_rate = 1 / denominator

# --- Step 4: Print the full calculation and the final answer ---
print("The firing rate (r) is calculated using the formula:")
print("r = 1 / (τ_ref + τ * ln((μ - V_reset) / (μ - V_th)))")
print("\nSubstituting the given values into the equation:")

print(f"μ (mean potential) = {V_ext * 1000} mV")
print(f"V_th (threshold) = {V_th * 1000} mV")
print(f"V_reset (reset) = {V_reset * 1000} mV")
print(f"τ (time constant) = {tau * 1000} ms")
print(f"τ_ref (refractory period) = {tau_ref * 1000} ms")

print("\nStep-by-step calculation:")
print(f"r = 1 / ({tau_ref} s + {tau} s * ln(({V_ext} V - {V_reset} V) / ({V_ext} V - {V_th} V)))")
print(f"r = 1 / ({tau_ref} s + {tau} s * ln({log_argument}))")
print(f"r = 1 / ({tau_ref} s + {tau} s * {natural_log:.6f})")
print(f"r = 1 / ({denominator:.6f} s)")
print(f"\nCalculated firing rate = {firing_rate:.2f} Hz")

final_answer = round(firing_rate)
print(f"\nThe firing rate of a typical neuron, rounded to the nearest integer, is: {final_answer} Hz")

<<<63>>>
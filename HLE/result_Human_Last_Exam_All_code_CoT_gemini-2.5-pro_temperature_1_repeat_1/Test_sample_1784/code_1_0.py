import math

# Step 1: Define parameters from the problem description.
# All voltage and time parameters are defined in mV and ms for clarity in the printout.
tau_m_ms = 20.0      # Membrane time constant in ms
V_reset_mv = 10.0    # Voltage reset in mV
V_th_mv = 20.0       # Voltage threshold in mV
tau_ref_ms = 2.0     # Refractory period in ms
mu_ext_mv = 30.0     # External input in mV

# Network parameters used to determine the mean input.
K_E = 1000
J_E_mv = 0.1
g = 4
K_I = 0.25 * K_E
J_I_mv = g * J_E_mv

# Step 2: Calculate the total mean input voltage (μ).
# The mean input from the network is proportional to (K_E * J_E - K_I * J_I).
# (1000 * 0.1) - (250 * 0.4) = 100 - 100 = 0.
# The network is balanced, so the mean recurrent input is zero.
# The total mean input μ is therefore equal to the external input.
mu_mv = mu_ext_mv

# Step 3: Use the appropriate formula for the supra-threshold regime.
# Since mu_mv (30 mV) > V_th_mv (20 mV), we use the deterministic formula.
# ν = 1 / (τ_ref + τ_m * log((μ - V_reset) / (μ - V_th)))

# Convert time units to seconds for calculation in Hz.
tau_m_s = tau_m_ms / 1000.0
tau_ref_s = tau_ref_ms / 1000.0

# Step 4: Calculate the final firing rate.
log_argument = (mu_mv - V_reset_mv) / (mu_mv - V_th_mv)
firing_rate_hz = 1.0 / (tau_ref_s + tau_m_s * math.log(log_argument))

# Step 5: Print the explanation, equation, and result as requested.
print("The firing rate (ν) of the neuron is calculated based on its properties and the total mean input voltage (μ).")
print(f"The total mean input voltage μ is {mu_mv} mV, which is above the threshold V_th of {V_th_mv} mV.")
print("\nUsing the formula for the supra-threshold firing rate:")
print("ν = 1 / (τ_ref + τ_m * log((μ - V_reset) / (μ - V_th)))")

print("\nPlugging in the given numbers:")
print(f"ν = 1 / ({tau_ref_s} s + {tau_m_s} s * log(({mu_mv} mV - {V_reset_mv} mV) / ({mu_mv} mV - {V_th_mv} mV)))")
print(f"ν = 1 / ({tau_ref_s} + {tau_m_s} * log({log_argument}))")

print("\nThe final firing rate is (in Hz):")
print(int(round(firing_rate_hz)))
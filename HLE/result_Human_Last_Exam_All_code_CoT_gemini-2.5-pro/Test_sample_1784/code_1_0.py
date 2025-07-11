import math

# Step 1: Define parameters in consistent units (ms, mV)
# Membrane time constant in ms
tau_m = 20.0
# Excitatory synaptic efficacy in mV
J_E = 0.1
# Voltage reset in mV
V_reset = 10.0
# Voltage threshold in mV
V_th = 20.0
# Refractory period in ms
tau_ref = 2.0
# Relative strength of inhibition to excitation
g = 4.0
# Number of excitatory connections per neuron
K_E = 1000.0
# Proportion of inhibitory connections relative to excitatory
p_I = 0.25
# External input in mV
V_ext = 30.0

# Step 2: Calculate derived parameters for the inhibitory population
# Number of inhibitory connections
K_I = p_I * K_E
# Inhibitory synaptic efficacy in mV
J_I = g * J_E

# Step 3: Analyze the network balance
# In a recurrent network, the mean potential drive (mu_V) is given by:
# mu_V = V_ext + tau_m * nu * (K_E * J_E - K_I * J_I)
# where 'nu' is the network firing rate. Let's check the balance term.
recurrent_excitation_term = K_E * J_E
recurrent_inhibition_term = K_I * J_I
net_recurrent_term = recurrent_excitation_term - recurrent_inhibition_term

# Because the net recurrent term is zero, the network is "balanced".
# This simplifies the mean potential drive to be equal to the external input.
mu_V = V_ext

print("Analyzing the network balance:")
print(f"Total recurrent excitation strength (K_E * J_E): {K_E} * {J_E} = {recurrent_excitation_term:.1f} mV")
print(f"Total recurrent inhibition strength (K_I * J_I): ({p_I} * {K_E}) * ({g} * {J_E}) = {K_I} * {J_I} = {recurrent_inhibition_term:.1f} mV")
print(f"Net recurrent input strength: {recurrent_excitation_term:.1f} - {recurrent_inhibition_term:.1f} = {net_recurrent_term:.1f} mV")
print("\nThe network is perfectly balanced, so the mean input potential is determined by the external input.")
print(f"\nMean Potential Drive (μ_V) = External Input (V_ext) = {mu_V:.1f} mV")

# Step 4: Calculate the firing rate
# Compare mean drive to threshold to confirm the firing regime.
print(f"Voltage Threshold (V_th) = {V_th:.1f} mV")
if mu_V > V_th:
    print("Since μ_V > V_th, the neuron is in the supra-threshold (mean-driven) firing regime.")

    # The formula for the firing rate (nu) in this regime is:
    # nu = 1 / (tau_ref + tau_m * ln((mu_V - V_reset) / (mu_V - V_th)))
    # We multiply by 1000 to convert from spikes/ms to spikes/s (Hz).
    
    numerator = mu_V - V_reset
    denominator = mu_V - V_th
    log_term = math.log(numerator / denominator)
    time_to_threshold = tau_m * log_term
    
    # Inter-spike interval in ms
    ISI_ms = tau_ref + time_to_threshold
    
    # Firing rate in Hz
    firing_rate_hz = 1000.0 / ISI_ms
    
    print("\nCalculating the firing rate (ν):")
    print("ν = 1000 / (τ_ref + τ_m * ln((μ_V - V_reset) / (μ_V - V_th)))")
    print("Plugging in the numbers:")
    print(f"ν = 1000 / ({tau_ref:.1f} + {tau_m:.1f} * ln(({mu_V:.1f} - {V_reset:.1f}) / ({mu_V:.1f} - {V_th:.1f})))")
    print(f"ν = 1000 / ({tau_ref:.1f} + {tau_m:.1f} * ln({numerator:.1f} / {denominator:.1f}))")
    print(f"ν = 1000 / ({tau_ref:.1f} + {tau_m:.1f} * {log_term:.4f})")
    print(f"ν = 1000 / ({tau_ref:.1f} + {time_to_threshold:.4f})")
    print(f"ν = 1000 / {ISI_ms:.4f}")
    print(f"ν ≈ {firing_rate_hz:.2f} Hz")
    
    # Final answer as a rounded integer
    final_answer = int(round(firing_rate_hz))
    print(f"\nThe final firing rate, rounded to the nearest integer, is: {final_answer} Hz.")

else:
    # This case is not expected based on the problem's parameters
    print("The neuron is in the sub-threshold regime, and its firing rate would be zero without fluctuations.")
    final_answer = 0

print(f"\n<<<{int(round(firing_rate_hz))}>>>")
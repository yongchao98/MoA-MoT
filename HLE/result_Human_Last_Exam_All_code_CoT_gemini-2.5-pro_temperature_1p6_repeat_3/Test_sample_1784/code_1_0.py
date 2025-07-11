import math

# Step 1: Define all parameters from the problem statement in consistent units (ms, mV)
tau_ms = 20.0       # Membrane time constant (ms)
J_E_mV = 0.1        # Excitatory synaptic efficacy (mV)
V_reset_mV = 10.0   # Voltage reset (mV)
V_th_mV = 20.0      # Voltage threshold (mV)
t_ref_ms = 2.0      # Refractory period (ms)
g = 4.0             # Relative inhibition to excitation
K_E = 1000          # Number of excitatory connections
inh_proportion = 0.25 # Proportion of inhibitory to excitatory connections
V_ext_mV = 30.0     # External input (mV)

# Step 2: Analyze the network inputs
# Calculate the number of inhibitory connections
K_I = K_E * inh_proportion

# The term "relative inhibition to excitation is 4" (g=4) implies that the inhibitory
# synaptic efficacy J_I is g times the excitatory one J_E.
J_I_mV = g * J_E_mV

# The mean input from the recurrent network connections depends on the firing rate 'r'
# and the balance between excitation and inhibition.
# Total recurrent excitatory drive strength: K_E * J_E
recurrent_excitation = K_E * J_E_mV
# Total recurrent inhibitory drive strength: K_I * J_I
recurrent_inhibition = K_I * J_I_mV

# If recurrent_excitation == recurrent_inhibition, the network is balanced.
# Let's check:
# recurrent_inhibition = (1000 * 0.25) * (4 * 0.1) = 250 * 0.4 = 100 mV
# recurrent_excitation = 1000 * 0.1 = 100 mV
# They are equal, so the average recurrent input is zero.

# The total mean voltage drive (mu_V) is therefore equal to the external input.
mu_V_mV = V_ext_mV

# Step 3: Calculate the firing rate
# We observe that the mean input mu_V (30 mV) is greater than the threshold V_th (20 mV).
# This means the neuron is in a mean-driven regime.
# We can approximate the firing rate by neglecting input noise and calculating the time
# it takes for the membrane potential to deterministically drift from V_reset to V_th.
# The formula for this drift time is: T_drift = tau * ln((mu_V - V_reset) / (mu_V - V_th))
# The total inter-spike interval (ISI) is T_drift + t_ref.
# The firing rate 'r' is 1 / ISI.

# Perform the calculation
if (mu_V_mV > V_th_mV):
    # The term inside the logarithm
    log_term = (mu_V_mV - V_reset_mV) / (mu_V_mV - V_th_mV)
    
    # Drift time in ms
    T_drift_ms = tau_ms * math.log(log_term)
    
    # Inter-spike interval in ms
    ISI_ms = t_ref_ms + T_drift_ms
    
    # Firing rate in Hz (1/s), so we convert ISI from ms to s by dividing by 1000.
    rate_Hz = 1000.0 / ISI_ms
    
    # Print the equation with values plugged in, as requested.
    # We use integers for values that are naturally so, and floats for the rest.
    print("The firing rate 'r' is calculated using the formula for a mean-driven integrate-and-fire neuron:")
    print(f"r = 1000 / (t_ref + τ * ln((μ_V - V_reset) / (μ_V - V_th)))")
    print("Plugging in the numbers (time in ms, voltage in mV):")
    print(f"r = 1000 / ({t_ref_ms:.0f} + {tau_ms:.0f} * ln(({mu_V_mV:.0f} - {V_reset_mV:.0f}) / ({mu_V_mV:.0f} - {V_th_mV:.0f})))")
    print(f"r = 1000 / ({t_ref_ms:.0f} + {tau_ms:.0f} * ln({log_term:.1f}))")
    print(f"r ≈ {rate_Hz:.2f} Hz")
    
    # Final answer as an integer
    final_rate_int = int(round(rate_Hz))
    print("\nThe final firing rate as an integer is:")
    print(final_rate_int)

else:
    # This case does not apply here but is included for completeness
    final_rate_int = 0
    print("Mean input is not above threshold; the firing rate is 0 Hz in this simplified model.")

<<<63>>>
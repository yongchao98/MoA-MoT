import math

# 1. Define the parameters from the problem description
tau_m_ms = 20.0  # Membrane time constant in ms
J_mv = 0.1     # Synaptic efficacy in mV
V_reset_mv = 10.0 # Voltage reset in mV
V_th_mv = 20.0   # Voltage threshold in mV
tau_ref_ms = 2.0 # Refractory period in ms
g = 4.0          # Relative inhibition to excitation
K_E = 1000.0     # Number of excitatory connections
inh_proportion = 0.25 # Proportion of inhibitory connections relative to excitatory
V_ext_mv = 30.0  # External input in mV

# 2. Calculate the number of inhibitory connections
K_I = K_E * inh_proportion

# 3. Check for balance between excitation and inhibition from recurrent connections.
# The mean input from the network is proportional to the firing rate ν and the term (K_E - g * K_I).
balance_factor = K_E - g * K_I

print("Step 1: Analyzing the network's recurrent connections.")
print(f"Number of excitatory connections (K_E): {int(K_E)}")
print(f"Number of inhibitory connections (K_I = K_E * {inh_proportion}): {int(K_I)}")
print(f"Relative inhibitory strength (g): {int(g)}")
print(f"The balance factor is K_E - g * K_I = {int(K_E)} - {int(g)} * {int(K_I)} = {balance_factor}")
print("\nSince the balance factor is 0, the network's recurrent connections are perfectly balanced.")
print("This means the average excitatory and inhibitory inputs cancel each other out.")

# 4. Determine the total mean input potential (mu).
# Because the recurrent input is balanced, the mean input potential (mu) is simply the external input.
mu_mv = V_ext_mv
print(f"\nStep 2: Determining the mean input potential to a neuron.")
print(f"The mean input potential (μ) is equal to the external input: {int(mu_mv)} mV.")

# 5. Calculate the firing rate using the standard formula.
# We can use units of ms for time and mV for voltage consistently.
# The formula for firing rate (nu) in Hz is: 1000 / (tau_ref_ms + tau_m_ms * log((mu_mv - V_reset_mv) / (mu_mv - V_th_mv)))
print(f"\nStep 3: Calculating the firing rate (ν).")
print("The formula for the firing rate of an integrate-and-fire neuron is:")
print("ν = 1 / (τ_ref + τ * ln((μ - V_reset) / (μ - V_th)))")

# Using values in ms and mV, we multiply the result by 1000 to get Hz
log_argument = (mu_mv - V_reset_mv) / (mu_mv - V_th_mv)
log_result = math.log(log_argument)
denominator = tau_ref_ms + tau_m_ms * log_result
firing_rate_hz = 1000.0 / denominator

# Display the final equation with all numbers plugged in
equation_str = (
    f"ν = 1000 / ({tau_ref_ms} + {tau_m_ms} * ln(({mu_mv} - {V_reset_mv}) / ({mu_mv} - {V_th_mv})))"
)
print("\nPlugging in the values:")
print(equation_str)

# Final result
final_rate = int(round(firing_rate_hz))
print(f"\nThe calculated firing rate is approximately {firing_rate_hz:.2f} Hz.")
print(f"The final answer, rounded to the nearest integer, is: {final_rate}")

<<<63>>>
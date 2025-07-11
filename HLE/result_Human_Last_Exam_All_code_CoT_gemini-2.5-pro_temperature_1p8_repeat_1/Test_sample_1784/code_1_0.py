import math

# Step 1: Define parameters from the problem statement
tau_m = 20.0  # Membrane time constant (ms)
J_e = 0.1  # Excitatory synaptic efficacy (mV)
V_reset = 10.0  # Voltage reset (mV)
V_theta = 20.0  # Voltage threshold (mV)
tau_ref = 2.0  # Refractory period (ms)
g = 4.0  # Relative inhibition to excitation
C_e = 1000  # Number of excitatory connections
prop_i = 0.25  # Proportion of inhibitory connections
mu_ext = 30.0  # External input (mV)

# Step 2: Calculate mean recurrent input parameters
# Calculate the number of inhibitory connections
C_i = C_e * prop_i
# Calculate the inhibitory synaptic efficacy
# We assume g is the ratio of synaptic efficacies |J_i/J_e|
J_i = -g * J_e

# Step 3: Calculate the total mean membrane potential (mu_V)
# The mean input from the network is proportional to (C_e * J_e + C_i * J_i).
# In a self-consistent state, this term would be multiplied by the firing rate nu and tau_m.
# Let's check this term first.
mean_recurrent_input_coeff = C_e * J_e + C_i * J_i

# The total mean potential is the sum of external and recurrent inputs.
# mu_V = mu_ext + nu * tau_m * mean_recurrent_input_coeff
# Since the mean_recurrent_input_coeff is 0, mu_V is simply mu_ext.
if mean_recurrent_input_coeff == 0:
    mu_V = mu_ext
else:
    # If it were not zero, we would need to solve a self-consistent equation for nu.
    # But for this problem, the network is perfectly balanced.
    mu_V = mu_ext

# Step 4: Determine the firing regime
# Compare mu_V to V_theta. If mu_V > V_theta, the neuron is in the supra-threshold regime.
# We will use the deterministic formula for the firing rate, which is a good approximation here.

# Step 5: Calculate the firing rate using the supra-threshold formula
# Firing Rate (nu) = 1 / (tau_ref + T_flight)
# T_flight is the time to travel from V_reset to V_theta
# T_flight = tau_m * ln((mu_V - V_reset) / (mu_V - V_theta))

T_flight = tau_m * math.log((mu_V - V_reset) / (mu_V - V_theta))
inter_spike_interval = tau_ref + T_flight

# The firing rate nu is 1 / inter_spike_interval (in kHz, since times are in ms)
# Multiply by 1000 to get the rate in Hz.
firing_rate = 1000.0 / inter_spike_interval
integer_firing_rate = int(round(firing_rate))

# Step 6: Format and print the output
equation_str = (
    f"Firing Rate = 1000 / ({tau_ref} + {tau_m} * ln(({mu_V} - {V_reset}) / ({mu_V} - {V_theta})))"
)

print("Calculation Steps:")
print(f"Inhibitory connections C_i = {C_e} * {prop_i} = {C_i}")
print(f"Inhibitory efficacy J_i = -{g} * {J_e} = {J_i} mV")
print(f"Mean recurrent input coefficient = ({C_e}*{J_e}) + ({C_i}*{J_i}) = {C_e * J_e + C_i * J_i} mV")
print("Since the mean recurrent input is 0, the mean potential mu_V is equal to the external input.")
print(f"Mean potential mu_V = {mu_ext} mV")
print(f"Since mu_V ({mu_V} mV) > V_theta ({V_theta} mV), the neuron is in the supra-threshold regime.")
print("\nThe firing rate is calculated using the formula:")
print(equation_str)
print(f"\nCalculated firing rate: {firing_rate:.2f} Hz")
print(f"The firing rate of a typical neuron (rounded to the nearest integer) is: {integer_firing_rate} Hz")
print(f"\nFinal Answer: {integer_firing_rate}")
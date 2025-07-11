import math

# Step 1: Define all given parameters in consistent units (mV, ms, etc.)
tau_m = 20.0  # Membrane time constant (ms)
J_E = 0.1     # Excitatory synaptic efficacy (mV)
V_reset = 10.0  # Voltage reset (mV)
V_thresh = 20.0 # Voltage threshold (mV)
tau_ref = 2.0   # Refractory period (ms)
g = 4.0       # Relative inhibition to excitation strength
K_E = 1000    # Number of excitatory connections
p_I = 0.25    # Proportion of inhibitory connections relative to excitatory
mu_ext = 30.0   # External input drive (mV)

# Step 2: Calculate parameters for the inhibitory population
K_I = K_E * p_I  # Number of inhibitory connections
J_I = J_E * g    # Inhibitory synaptic efficacy (mV)

# Step 3: Determine the mean input voltage (μ)
# The total drive from excitatory connections is K_E * J_E
# The total drive from inhibitory connections is K_I * J_I
drive_E = K_E * J_E
drive_I = K_I * J_I

print(f"Number of inhibitory connections (K_I): {int(K_I)}")
print(f"Inhibitory synaptic efficacy (J_I): {J_I:.1f} mV")
print(f"Total excitatory drive per spike: {drive_E:.1f}")
print(f"Total inhibitory drive per spike: {drive_I:.1f}")
print("-" * 20)

# Because the excitatory and inhibitory drives are equal, the mean recurrent
# network input is zero. The net mean input is therefore just the external input.
mu = mu_ext
print(f"The network is balanced. Mean input voltage (μ) is equal to the external input: {mu:.1f} mV")
print("-" * 20)

# Step 4: Calculate the firing rate (ν) using the LIF neuron formula
# The formula for firing rate is 1 / (tau_ref + T), where T is the time to reach threshold.
# T = tau_m * ln((μ - V_reset) / (μ - V_thresh))
# We are calculating the rate in Hz, so we convert time from ms to s by dividing by 1000.
# Or, calculate rate in spikes/ms and multiply by 1000. Let's do that.

# Calculate time to charge from reset to threshold
time_to_charge = tau_m * math.log((mu - V_reset) / (mu - V_thresh))

# Total inter-spike interval (ISI) in ms
isi = tau_ref + time_to_charge

# Firing rate in Hz (1000 ms / 1 s)
fire_rate = 1000 / isi

# Print the final calculation step-by-step
print("Calculating the firing rate:")
print(f"Time to charge from V_reset to V_thresh (T) = {tau_m:.0f} * ln(({mu:.0f} - {V_reset:.0f}) / ({mu:.0f} - {V_thresh:.0f}))")
print(f"T = {tau_m:.0f} * ln({mu - V_reset:.0f} / {mu - V_thresh:.0f}) = {time_to_charge:.2f} ms")
print(f"Inter-spike Interval (ISI) = τ_ref + T = {tau_ref:.0f} ms + {time_to_charge:.2f} ms = {isi:.2f} ms")
print(f"Firing Rate (ν) = 1000 / ISI = 1000 / {isi:.2f} ms = {fire_rate:.2f} Hz")
print("-" * 20)

# Final answer rounded to the nearest integer
final_answer = round(fire_rate)
print(f"The firing rate of a typical neuron is {final_answer} Hz.")

# The problem requests a specific output format.
# Let's print the final numerical equation clearly.
print("\nFinal equation for the firing rate (ν) in Hz:")
print(f"ν = 1000 / ({tau_ref:.0f} + {tau_m:.0f} * ln(({mu:.0f} - {V_reset:.0f}) / ({mu:.0f} - {V_thresh:.0f})))")

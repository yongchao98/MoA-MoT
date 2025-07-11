import math

def calculate_neuron_firing_rate():
    """
    Calculates the firing rate of a neuron in a balanced random network.
    """
    # 1. Define the given parameters
    tau_m = 20.0  # Membrane time constant (ms)
    J_E = 0.1  # Excitatory synaptic efficacy (mV)
    V_reset = 10.0  # Voltage reset (mV)
    V_th = 20.0  # Voltage threshold (mV)
    tau_ref = 2.0  # Refractory period (ms)
    g = 4.0  # Relative inhibition to excitation
    K_E = 1000.0  # Number of excitatory connections
    K_I_prop = 0.25  # Proportion of inhibitory to excitatory connections
    V_ext = 30.0  # External input driving potential (mV)

    # 2. Analyze the network input
    # Calculate the number of inhibitory connections and their strength
    K_I = K_E * K_I_prop
    J_I = g * J_E
    
    # The mean recurrent input is proportional to (K_E * J_E - K_I * J_I).
    # (1000 * 0.1) - (250 * 0.4) = 100 - 100 = 0.
    # The network is balanced, so the mean recurrent input is zero.
    # Therefore, the total mean input potential 'mu' is the external input.
    mu = V_ext

    # 3. Apply the firing rate formula if input is above threshold
    if mu <= V_th:
        firing_rate_hz = 0
    else:
        # Time for membrane potential to go from V_reset to V_th (in ms)
        time_to_threshold = tau_m * math.log((mu - V_reset) / (mu - V_th))
        
        # Total time between spikes is the sum of the time to threshold and refractory period (in ms)
        total_period_ms = tau_ref + time_to_threshold
        
        # Firing rate in Hz (1000 ms / 1 s)
        firing_rate_hz = 1000.0 / total_period_ms
    
    # 4. Print the equation and the final integer answer
    print("The firing rate ν (in Hz) is calculated using the formula:")
    print("ν = 1000 / (τ_ref + τ_m * ln((μ - V_reset) / (μ - V_th)))")
    print("\nPlugging in the numbers:")
    # Print the equation with each number explicitly shown
    print(f"ν = 1000 / ({tau_ref} + {tau_m} * ln(({mu} - {V_reset}) / ({mu} - {V_th})))")
    
    final_answer = int(round(firing_rate_hz))
    print(f"\nCalculation result: {firing_rate_hz:.2f} Hz")
    print(f"The firing rate of a typical neuron is {final_answer} Hz.")

# Execute the function to get the answer
calculate_neuron_firing_rate()
<<<63>>>
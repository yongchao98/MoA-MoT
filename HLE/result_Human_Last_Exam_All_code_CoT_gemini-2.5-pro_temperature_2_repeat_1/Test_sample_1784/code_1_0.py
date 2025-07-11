import math

def calculate_firing_rate():
    """
    Calculates the firing rate of a neuron in a balanced random network based on given parameters.
    """
    # 1. Define Parameters in SI units (Volts, Seconds)
    tau_ms = 20.0        # Membrane time constant (ms)
    J_mv = 0.1           # Synaptic efficacy (mV)
    V_reset_mv = 10.0    # Reset potential (mV)
    V_thresh_mv = 20.0   # Threshold potential (mV)
    tau_ref_ms = 2.0     # Refractory period (ms)
    g = 4.0              # Relative inhibition strength
    K_E = 1000           # Number of excitatory connections
    K_I_ratio = 0.25     # Proportion of inhibitory connections
    V_ext_mv = 30.0      # External input (mV)

    # 2. Calculate the properties of the recurrent network input
    K_I = K_E * K_I_ratio
    
    # In mean-field theory, the mean recurrent input potential is proportional to:
    # (K_E * J * ν * τ) - (g * K_I * J * ν * τ) = (K_E - g * K_I) * J * ν * τ
    # Let's check the balance factor
    balance_factor = K_E - g * K_I
    
    print("--- Firing Rate Calculation ---")
    print("\nStep 1: Determine the mean membrane potential (μ_V).")
    print(f"The network's excitatory drive is proportional to K_E = {K_E}.")
    print(f"The effective inhibitory drive is proportional to g * K_I = {g} * {K_I} = {g * K_I}.")
    
    if balance_factor == 0:
        print("The recurrent excitatory and inhibitory inputs are perfectly balanced and cancel each other out.")
        print(f"Therefore, the mean membrane potential μ_V is equal to the external input V_ext.")
        mu_V_mv = V_ext_mv
        print(f"μ_V = {mu_V_mv:.0f} mV\n")
    else:
        # This case is not applicable for the problem but included for completeness.
        print("The network is not perfectly balanced. This would require solving a self-consistent equation.")
        # For this problem, we proceed with the balanced case.
        mu_V_mv = V_ext_mv

    # 3. Calculate the firing rate using the LIF neuron firing rate formula.
    # We check if the mean potential is above threshold.
    if mu_V_mv <= V_thresh_mv:
        nu_hz = 0
        print("Mean potential is not above threshold; the firing rate is 0 Hz.")
    else:
        print("Step 2: Substitute values into the firing rate equation.")
        print("The time between spikes is given by T = τ_ref + τ * log((μ_V - V_reset) / (μ_V - V_thresh)).")
        print("The firing rate is ν = 1000 / T (where T is in ms).\n")
        
        print("The final equation with values (in ms and mV):")
        # To avoid division by zero if mu_V_mv == V_thresh_mv
        log_argument = (mu_V_mv - V_reset_mv) / (mu_V_mv - V_thresh_mv)
        
        print(f"ν (Hz) = 1000 / ({tau_ref_ms:.0f} + {tau_ms:.0f} * log(({mu_V_mv:.0f} - {V_reset_mv:.0f}) / ({mu_V_mv:.0f} - {V_thresh_mv:.0f})))")
        print(f"ν (Hz) = 1000 / ({tau_ref_ms:.0f} + {tau_ms:.0f} * log({mu_V_mv - V_reset_mv:.0f} / {mu_V_mv - V_thresh_mv:.0f}))")
        print(f"ν (Hz) = 1000 / ({tau_ref_ms:.0f} + {tau_ms:.0f} * log({log_argument}))")
        
        # Perform the actual calculation
        total_period_ms = tau_ref_ms + tau_ms * math.log(log_argument)
        nu_hz = 1000 / total_period_ms
        
        print(f"ν (Hz) = 1000 / ({total_period_ms:.2f})")
        print(f"ν (Hz) = {nu_hz:.2f}\n")
        
    # 4. Print the final answer rounded to the nearest integer.
    final_rate = int(round(nu_hz))
    print(f"The firing rate rounded to the nearest integer is: {final_rate}")

# Execute the function
calculate_firing_rate()
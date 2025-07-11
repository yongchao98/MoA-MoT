import math

def calculate_firing_rate():
    """
    Calculates the firing rate of a neuron in a balanced recurrent network.
    """
    # 1. Define the parameters from the problem description (units: ms, mV)
    tau_m = 20.0        # Membrane time constant (ms)
    V_reset = 10.0      # Reset voltage (mV)
    V_th = 20.0         # Voltage threshold (mV)
    tau_ref = 2.0       # Refractory period (ms)
    g = 4.0             # Relative inhibition to excitation
    gamma = 0.25        # Ratio of inhibitory to excitatory connections
    V_ext = 30.0        # External input (mV)

    # 2. Check for network balance and determine mean input potential 'mu'
    # The mean input 'mu' is the sum of recurrent and external inputs.
    # The recurrent part is proportional to (1 - g * gamma).
    balance_check = g * gamma
    # Since g * gamma = 4.0 * 0.25 = 1.0, the recurrent inputs cancel out.
    # Therefore, the mean potential is simply the external input.
    mu = V_ext

    # 3. Calculate firing rate using the supra-threshold formula for a LIF neuron.
    # Firing rate nu = 1 / T_ISI, where T_ISI is the inter-spike interval.
    # The rate is converted to Hz by using 1000 in the numerator (since time is in ms).
    # nu(Hz) = 1000 / (tau_ref + tau_m * ln((mu - V_reset) / (mu - V_th)))
    
    # Calculate the components of the formula for clarity
    term_mu_vreset = mu - V_reset
    term_mu_vth = mu - V_th
    log_arg = term_mu_vreset / term_mu_vth
    
    # Calculate the Inter-Spike Interval (ISI) in milliseconds
    T_ISI_ms = tau_ref + tau_m * math.log(log_arg)

    # Calculate the firing rate in Hz
    firing_rate_hz = 1000 / T_ISI_ms
    final_answer_int = round(firing_rate_hz)

    # 4. Print the calculation step-by-step
    print("The network balance factor g * γ is calculated as:")
    print(f"g * γ = {g} * {gamma} = {balance_check}")
    print("Since the result is 1, the network is in a balanced state. The average recurrent inputs cancel out, so the mean input potential μ is equal to the external input V_ext.")
    print(f"μ = {V_ext} mV\n")
    print("Since μ > V_th (30 mV > 20 mV), we use the supra-threshold LIF neuron firing rate formula:")
    print("ν = 1000 / (τ_ref + τ_m * ln((μ - V_reset) / (μ - V_th)))\n")
    print("Plugging in the numbers:")
    print(f"ν = 1000 / ({tau_ref} + {tau_m} * ln(({mu} - {V_reset}) / ({mu} - {V_th})))")
    print(f"ν = 1000 / ({tau_ref} + {tau_m} * ln({term_mu_vreset} / {term_mu_vth}))")
    print(f"ν = 1000 / ({tau_ref} + {tau_m} * ln({log_arg}))")
    print(f"ν = 1000 / ({T_ISI_ms:.2f})")
    print(f"ν ≈ {firing_rate_hz:.2f} Hz\n")
    print("The final firing rate, rounded to the nearest integer, is:")
    print(f"{final_answer_int}")

# Execute the function to get the answer
calculate_firing_rate()
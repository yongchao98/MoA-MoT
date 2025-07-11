import math

def calculate_firing_rate():
    """
    Calculates the firing rate of a neuron in a balanced random network.
    """
    # 1. Define Parameters in consistent units (ms and mV)
    tau = 20.0  # Membrane time constant (ms)
    J = 0.1     # Excitatory synaptic efficacy (mV)
    V_reset = 10.0 # Voltage reset (mV)
    V_th = 20.0    # Voltage threshold (mV)
    tau_ref = 2.0  # Refractory period (ms)
    g = 4.0        # Relative inhibition to excitation (J_I / J_E)
    K_E = 1000     # Number of excitatory connections
    K_I = 250      # Number of inhibitory connections (0.25 * K_E)
    V_ext = 30.0   # External input (mV)

    # 2. Calculate Mean Membrane Potential (mu_V)
    # The network is balanced: K_E * J = 1000 * 0.1 = 100 mV
    # and K_I * (g*J) = 250 * (4*0.1) = 250 * 0.4 = 100 mV.
    # The mean recurrent input is proportional to (K_E*J - K_I*J_I) which is 0.
    # Therefore, the mean potential is determined by the external input.
    mu_V = V_ext

    print(f"Neuron and Network Parameters:")
    print(f"Membrane time constant (τ): {tau} ms")
    print(f"Voltage threshold (V_th): {V_th} mV")
    print(f"Voltage reset (V_reset): {V_reset} mV")
    print(f"Refractory period (τ_ref): {tau_ref} ms")
    print(f"External input leads to a mean potential (μ_V) of: {mu_V} mV")
    print("-" * 30)
    
    # 3. Determine Firing Regime
    # Since mu_V (30 mV) > V_th (20 mV), the neuron is in the supra-threshold regime.

    # 4. Calculate Inter-Spike Interval (T_isi)
    # The formula is: T_isi = τ * ln((μ_V - V_reset) / (μ_V - V_th))
    # Check if the argument of the logarithm is positive
    if (mu_V - V_th) <= 0:
        print("Neuron is not in the supra-threshold regime. This formula cannot be applied.")
        # In a sub-threshold case, rate would be 0 without fluctuations.
        rate_hz = 0
    else:
        # Calculate argument of the logarithm
        log_arg = (mu_V - V_reset) / (mu_V - V_th)
        
        # Calculate T_isi
        T_isi = tau * math.log(log_arg)
        
        # 5. Calculate total period and firing rate
        total_period_ms = tau_ref + T_isi
        # Rate in Hz = 1 / (period in seconds) = 1000 / (period in ms)
        rate_hz = 1000 / total_period_ms
        
        # Print the detailed calculation steps
        print("Calculating the firing rate (ν):")
        print(f"ν = 1 / (τ_ref + τ * ln((μ_V - V_reset) / (μ_V - V_th)))")
        print(f"ν = 1 / ({tau_ref} ms + {tau} ms * ln(({mu_V} mV - {V_reset} mV) / ({mu_V} mV - {V_th} mV)))")
        print(f"ν = 1 / ({tau_ref} ms + {tau} ms * ln({log_arg:.2f}))")
        print(f"ν = 1 / ({tau_ref} ms + {T_isi:.2f} ms)")
        print(f"ν = 1 / ({total_period_ms:.2f} ms)")
        print(f"ν = {rate_hz:.2f} Hz")
        print("-" * 30)

    # 6. Final Answer as an integer
    final_rate_int = round(rate_hz)
    print(f"The firing rate of a typical neuron is {final_rate_int} Hz.")

if __name__ == '__main__':
    calculate_firing_rate()
    # To conform to the output format, we extract the final integer answer
    # This calculation is repeated here just to isolate the final numeric answer for the format requirement.
    tau, V_reset, V_th, tau_ref, mu_V = 20.0, 10.0, 20.0, 2.0, 30.0
    T_isi = tau * math.log((mu_V - V_reset) / (mu_V - V_th))
    total_period_ms = tau_ref + T_isi
    rate_hz = 1000 / total_period_ms
    final_rate_int = round(rate_hz)
    print(f'<<<{final_rate_int}>>>')
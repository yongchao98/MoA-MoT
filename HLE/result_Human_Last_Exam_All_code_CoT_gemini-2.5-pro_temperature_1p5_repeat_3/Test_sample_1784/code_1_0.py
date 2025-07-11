import math

def calculate_firing_rate():
    """
    Calculates the firing rate of a neuron in a randomly connected network
    based on the provided biophysical and network parameters.
    """
    # Step 1: Define all the given parameters in consistent units (ms and mV).
    tau = 20.0        # Membrane time constant (ms)
    J_e = 0.1         # Excitatory synaptic efficacy (mV)
    V_reset = 10.0    # Voltage reset (mV)
    V_th = 20.0       # Voltage threshold (mV)
    tau_ref = 2.0     # Refractory period (ms)
    g = 4.0           # Relative strength of inhibition to excitation
    K_e = 1000        # Number of excitatory connections per neuron
    inh_prop = 0.25   # Proportion of inhibitory connections
    V_ext = 30.0      # External input driving potential (mV)

    # Step 2: Analyze the mean network input.
    # The network is balanced if the mean excitatory and inhibitory inputs cancel out.
    # Mean input is proportional to (K_e * J_e) - (K_i * J_i)
    K_i = K_e * inh_prop
    J_i = J_e * g
    
    # Check for balance:
    # Total excitatory weight: K_e * J_e = 1000 * 0.1 = 100
    # Total inhibitory weight: K_i * J_i = (1000 * 0.25) * (0.1 * 4) = 250 * 0.4 = 100
    # Since they are equal, the mean recurrent input is zero.
    # Therefore, the steady-state voltage (V_ss) is determined by the external input.
    V_ss = V_ext

    # Step 3: Determine the firing regime.
    # Since V_ss (30 mV) is greater than V_th (20 mV), the neuron is in the
    # supra-threshold firing regime.

    # Step 4 & 5: Apply the supra-threshold firing rate formula and calculate the rate.
    # The inter-spike interval (T_isi) is the sum of the refractory period and the time
    # it takes for the voltage to travel from V_reset to V_th.
    # T_isi = tau_ref + tau * ln((V_ss - V_reset) / (V_ss - V_th))
    # Firing rate (in Hz) = 1000 / T_isi (since time constants are in ms)

    numerator_log = V_ss - V_reset
    denominator_log = V_ss - V_th
    log_argument = numerator_log / denominator_log
    
    T_isi = tau_ref + tau * math.log(log_argument)
    firing_rate = 1000 / T_isi
    
    # Print the explanation and the step-by-step calculation of the formula.
    print("The firing rate (ν) is calculated using the formula for a supra-threshold integrate-and-fire neuron.")
    print("The formula is: ν = 1000 / (τ_ref + τ * ln((V_ss - V_reset) / (V_ss - V_th)))")
    print(f"\nPlugging in the numbers:")
    print(f"ν = 1000 / ({tau_ref} + {tau} * ln(({V_ss} - {V_reset}) / ({V_ss} - {V_th})))")
    print(f"ν = 1000 / ({tau_ref} + {tau} * ln({numerator_log} / {denominator_log}))")
    print(f"ν = 1000 / ({tau_ref} + {tau} * ln({log_argument}))")
    print(f"ν = 1000 / {T_isi:.2f}")
    print(f"ν ≈ {firing_rate:.2f} Hz")
    
    final_rate = int(round(firing_rate))
    print(f"\nThe firing rate of a typical neuron is {final_rate} Hz.")
    
    return final_rate

if __name__ == '__main__':
    calculate_firing_rate()
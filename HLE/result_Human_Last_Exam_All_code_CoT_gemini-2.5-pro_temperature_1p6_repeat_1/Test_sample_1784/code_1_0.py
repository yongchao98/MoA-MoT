import math

def calculate_neuron_firing_rate():
    """
    Calculates the firing rate of a neuron in a randomly connected network
    based on the provided parameters.
    """
    # 1. Define parameters from the problem description in consistent units (ms, mV)
    tau_m = 20.0      # Membrane time constant (ms)
    J = 0.1           # Excitatory synaptic efficacy (mV)
    V_reset = 10.0    # Reset voltage (mV)
    V_th = 20.0       # Voltage threshold (mV)
    tau_ref = 2.0     # Absolute refractory period (ms)
    g = 4.0           # Relative strength of inhibition to excitation
    K_E = 1000.0      # Number of excitatory connections
    K_I_prop = 0.25   # Proportion of inhibitory to excitatory connections
    mu_ext = 30.0     # External input (mV)

    # 2. Determine derived parameters
    K_I = K_E * K_I_prop  # Number of inhibitory connections
    J_E = J               # Efficacy of excitatory synapses
    J_I = g * J           # Efficacy of inhibitory synapses

    # 3. Calculate the net recurrent synaptic weight
    # This shows that the network is in a balanced state where recurrent
    # excitation and inhibition cancel each other out on average.
    net_recurrent_weight = K_E * J_E - K_I * J_I

    # The total mean input potential 'mu' is the sum of external and recurrent inputs.
    # In a balanced network (net_recurrent_weight=0), mu simplifies to mu_ext.
    mu = mu_ext

    # 4. Check if the neuron will fire
    if mu <= V_th:
        print(f"The mean input potential ({mu} mV) is not above the threshold ({V_th} mV), so the neuron will not fire.")
        firing_rate_hz = 0
    else:
        # 5. Calculate the firing rate using the LIF neuron firing rate formula.
        # The result of the formula is the rate in spikes per ms (kHz).
        
        # Calculate the denominator term by term for clarity
        delta_V_num = mu - V_reset
        delta_V_den = mu - V_th
        log_term = math.log(delta_V_num / delta_V_den)
        isi_no_refractory = tau_m * log_term
        
        # Inter-spike interval in ms
        T_isi_ms = tau_ref + isi_no_refractory
        
        # Firing rate in spikes per ms (kHz)
        firing_rate_khz = 1 / T_isi_ms
        
        # Convert firing rate to spikes per second (Hz)
        firing_rate_hz = firing_rate_khz * 1000
        
        # Print the equation and the values being plugged in
        print("The firing rate (v) of a leaky integrate-and-fire neuron is given by:")
        print("v = 1 / (t_ref + t_m * ln((mu - V_reset) / (mu - V_th)))")
        print("\nPlugging in the given values:")
        
        equation_str = (
            f"v = 1 / ({tau_ref} + {tau_m} * ln(({mu} - {V_reset}) / ({mu} - {V_th}))) * 1000 Hz"
        )
        print(equation_str)
        
        simplified_equation_str = (
            f"v = 1 / ({tau_ref} + {tau_m} * ln({delta_V_num / delta_V_den})) * 1000 Hz"
        )
        print(simplified_equation_str)


    # 6. Round the final answer to the nearest integer
    final_rate = int(round(firing_rate_hz))
    
    print(f"\nCalculated firing rate: {firing_rate_hz:.2f} Hz")
    print(f"The firing rate of a typical neuron, rounded to the nearest integer, is: {final_rate} Hz")

    return final_rate

# Execute the function
final_answer = calculate_neuron_firing_rate()
print(f'<<<{final_answer}>>>')

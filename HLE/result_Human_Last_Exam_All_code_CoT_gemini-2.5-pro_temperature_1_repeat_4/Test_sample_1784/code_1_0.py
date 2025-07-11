import math

def calculate_firing_rate():
    """
    Calculates the firing rate of a neuron in a balanced network based on given parameters.
    """
    # Neuron and network parameters (using mV for potential and ms for time)
    tau = 20.0          # Membrane time constant (ms)
    J_E = 0.1           # Excitatory synaptic efficacy (mV)
    V_reset = 10.0      # Reset potential (mV)
    V_th = 20.0         # Threshold potential (mV)
    tau_ref = 2.0       # Refractory period (ms)
    g = 4.0             # Relative strength of inhibition (J_I / J_E)
    K_E = 1000          # Number of excitatory connections
    inhib_prop = 0.25   # Proportion of inhibitory to excitatory connections
    V_ext = 30.0        # External input that sets the mean potential (mV)

    # Step 1: Determine inhibitory parameters and check for balance
    J_I = g * J_E
    K_I = inhib_prop * K_E
    
    # The total drive from recurrent connections is proportional to (J_E*K_E - J_I*K_I)
    # Excitatory drive: 0.1 mV * 1000 = 100 mV
    # Inhibitory drive: (4 * 0.1 mV) * (0.25 * 1000) = 0.4 mV * 250 = 100 mV
    # Since they are equal, the mean recurrent input is zero. The network is balanced.

    # Step 2: Determine total mean input potential
    # In a balanced network, the total mean input is equal to the external input.
    mu_total = V_ext

    # Step 3: Calculate the firing rate using the LIF neuron formula
    # Check if the neuron is in the supra-threshold regime (mu_total > V_th)
    if mu_total > V_th:
        # The equation for the firing rate (in Hz) is:
        # ν = 1000 / (τ_ref + τ * ln((μ_total - V_reset) / (μ_total - V_th)))
        # The factor of 1000 converts the result from 1/ms to 1/s (Hz).
        
        log_argument_numerator = mu_total - V_reset
        log_argument_denominator = mu_total - V_th
        log_argument = log_argument_numerator / log_argument_denominator

        T_charge = tau * math.log(log_argument)
        ISI = tau_ref + T_charge
        firing_rate = 1000.0 / ISI

        # Step 4: Output the results
        print("The firing rate (ν) is calculated using the formula for a supra-threshold Leaky Integrate-and-Fire neuron.")
        print("Final Equation:")
        print(f"ν = 1000 / (τ_ref + τ * ln((μ_total - V_reset) / (μ_total - V_th)))")
        print("\nSubstituting the numerical values:")
        print(f"ν = 1000 / ({tau_ref} + {tau} * ln(({mu_total} - {V_reset}) / ({mu_total} - {V_th})))")
        print(f"ν = 1000 / ({tau_ref} + {tau} * ln({log_argument_numerator} / {log_argument_denominator}))")
        print(f"ν = 1000 / ({tau_ref} + {tau} * {math.log(log_argument):.4f})")
        print(f"ν = 1000 / ({ISI:.4f})")
        
        final_answer = int(round(firing_rate))
        print(f"\nCalculated Firing Rate: {firing_rate:.2f} Hz")
        print(f"The firing rate of a typical neuron is {final_answer} Hz.")
        
        return final_answer
    else:
        # If input is not above threshold, firing rate is 0 in this model
        print("Input potential is not above threshold. Firing rate is 0 Hz.")
        return 0

# Execute the function and capture the final integer answer for submission
final_rate = calculate_firing_rate()
print(f'<<<{final_rate}>>>')
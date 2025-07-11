import math

def calculate_firing_rate():
    """
    Calculates the firing rate of a typical neuron in a randomly connected network
    based on the given parameters.
    """
    # 1. Define parameters in consistent units (mV, s)
    tau_m = 0.020  # Membrane time constant (20 ms)
    J = 0.1  # Excitatory synaptic efficacy (mV)
    V_reset = 10.0  # Voltage reset (mV)
    V_th = 20.0  # Voltage threshold (mV)
    tau_ref = 0.002  # Refractory period (2 ms)
    g = 4.0  # Relative inhibition to excitation
    K_E = 1000.0  # Number of excitatory connections
    inhibitory_proportion = 0.25 # Proportion of inhibitory connections
    mu_ext = 30.0  # External input (mV)

    # 2. Calculate derived network properties
    K_I = K_E * inhibitory_proportion  # Number of inhibitory connections
    J_E = J # Excitatory synaptic efficacy
    J_I = g * J  # Inhibitory synaptic efficacy

    print("Step 1: Calculating the total recurrent synaptic drive factor.")
    # This factor is multiplied by tau_m * v to get the recurrent potential
    recurrent_drive_factor = (K_E * J_E) - (K_I * J_I)
    print(f"Number of inhibitory connections (K_I) = {K_I}")
    print(f"Inhibitory synaptic efficacy (J_I) = {J_I:.1f} mV")
    print(f"Recurrent drive factor (K_E*J_E - K_I*J_I) = ({K_E}*{J_E}) - ({K_I}*{J_I}) = {recurrent_drive_factor:.1f} mV")
    print("-" * 30)

    print("Step 2: Calculating the mean membrane potential (mu_V).")
    # In a self-consistent model, mu_V = tau_m * v * recurrent_drive_factor + mu_ext
    # Since the recurrent drive factor is 0, the equation simplifies.
    mu_V = mu_ext
    print("The network is in a balanced state where the mean excitatory and inhibitory inputs cancel out.")
    print(f"Therefore, the mean membrane potential is determined by the external input.")
    print(f"mu_V = {mu_V:.1f} mV")
    print("-" * 30)

    print("Step 3: Calculating the firing rate (v).")
    # Check if the neuron will fire
    if mu_V <= V_th:
        firing_rate = 0
        print("Mean potential is not above threshold, so the firing rate is 0 Hz.")
    else:
        # Use the LIF neuron firing rate formula for a constant input
        log_term = math.log((mu_V - V_reset) / (mu_V - V_th))
        firing_rate = 1 / (tau_ref + tau_m * log_term)
        
        print("The firing rate 'v' is calculated using the formula:")
        print("v = 1 / (tau_ref + tau_m * log((mu_V - V_reset) / (mu_V - V_th)))")
        print("Plugging in the numbers:")
        print(f"v = 1 / ({tau_ref} + {tau_m} * log(({mu_V} - {V_reset}) / ({mu_V} - {V_th})))")
        print(f"v = 1 / ({tau_ref} + {tau_m} * log({(mu_V - V_reset)} / {(mu_V - V_th)}))")
        print(f"v = {firing_rate:.2f} Hz")
    
    print("-" * 30)
    
    # Final answer as an integer
    final_answer = int(round(firing_rate))
    print(f"The firing rate of a typical neuron, rounded to the nearest integer, is {final_answer} Hz.")
    
    return final_answer

if __name__ == '__main__':
    final_rate = calculate_firing_rate()
    print(f'<<<{final_rate}>>>')
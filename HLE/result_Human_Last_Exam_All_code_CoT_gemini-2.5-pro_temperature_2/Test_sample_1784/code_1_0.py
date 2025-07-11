import math

def calculate_firing_rate():
    """
    Calculates the firing rate of a neuron based on the provided network parameters.
    """
    # Step 1: Define parameters in base units (seconds, volts)
    tau_m = 20e-3       # Membrane time constant (s)
    J = 0.1e-3          # Base synaptic efficacy (V)
    V_reset = 10e-3     # Voltage reset (V)
    V_thresh = 20e-3    # Voltage threshold (V)
    tau_ref = 2e-3      # Refractory period (s)
    g = 4.0             # Relative strength of inhibition to excitation
    K_E = 1000          # Number of excitatory connections
    K_I_ratio = 0.25    # Proportion of inhibitory connections
    mu_V_ext = 30e-3    # External input driving potential (V)

    # Step 2: Determine synaptic strengths for excitatory and inhibitory neurons
    J_E = J
    J_I = J * g
    K_I = K_E * K_I_ratio

    # Step 3: Calculate the mean membrane potential (mu_V)
    # The mean drive from recurrent connections is proportional to (K_E * J_E - K_I * J_I)
    recurrent_drive_factor = (K_E * J_E) - (K_I * J_I)

    # Since recurrent_drive_factor is 0.0, the network is 'balanced', and the
    # mean potential is determined solely by the external input.
    mu_V = mu_V_ext

    # Step 4 & 5: Analyze the regime and calculate the firing rate.
    # Since mu_V (30 mV) > V_thresh (20 mV), the neuron is in the suprathreshold regime.
    # We use the formula for a noiseless LIF neuron as a reasonable approximation.
    if mu_V <= V_thresh:
        firing_rate = 0
        print("Mean potential is not above threshold. Firing rate is 0 Hz in the noise-free limit.")
    else:
        # T_isi = tau_m * ln((mu_V - V_reset) / (mu_V - V_thresh))
        # firing_rate = 1 / (tau_ref + T_isi)
        
        log_argument_num = mu_V - V_reset
        log_argument_den = mu_V - V_thresh
        log_argument = log_argument_num / log_argument_den
        
        T_isi = tau_m * math.log(log_argument)
        firing_rate = 1 / (tau_ref + T_isi)

        # Step 6: Output the calculation details
        print("The firing rate (ν) is calculated using the formula for a noiseless LIF neuron in the suprathreshold regime:")
        print("ν = 1 / (τ_ref + τ_m * ln((μ_V - V_reset) / (μ_V - V_thresh)))\n")
        print("Plugging in the given values (in Volts and seconds):")
        print(f"μ_V = {mu_V:.3f}, V_thresh = {V_thresh:.3f}, V_reset = {V_reset:.3f}, τ_m = {tau_m:.3f}, τ_ref = {tau_ref:.3f}\n")
        
        # We must output each number in the final equation.
        print(f"ν = 1 / ({tau_ref} + {tau_m} * ln(({mu_V} - {V_reset}) / ({mu_V} - {V_thresh})))")
        print(f"ν = 1 / ({tau_ref} + {tau_m} * ln({log_argument_num:.2f} / {log_argument_den:.2f}))")
        print(f"ν = 1 / ({tau_ref} + {tau_m} * ln({log_argument}))")
        print(f"ν = 1 / ({tau_ref:.3f} + {T_isi:.5f})")
        print(f"ν = 1 / {tau_ref + T_isi:.5f}")
        print(f"ν ≈ {firing_rate:.2f} Hz\n")
    
    final_answer = int(round(firing_rate))
    print(f"The final calculated firing rate, rounded to the nearest integer, is {final_answer} Hz.")

    # Store final answer for submission format
    global answer
    answer = final_answer


if __name__ == '__main__':
    calculate_firing_rate()
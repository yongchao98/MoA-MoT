import math

def calculate_firing_rate():
    """
    Calculates the firing rate of a neuron in a balanced network based on given parameters.
    """
    # --- Parameters (using ms and mV for clarity in output) ---
    tau_m_ms = 20.0      # Membrane time constant (ms)
    J_mV = 0.1         # Synaptic efficacy (mV)
    V_reset_mV = 10.0    # Voltage reset (mV)
    V_th_mV = 20.0       # Voltage threshold (mV)
    tau_ref_ms = 2.0     # Refractory period (ms)
    g = 4.0            # Relative inhibition to excitation
    K_E = 1000         # Number of excitatory connections
    prop_I = 0.25      # Proportion of inhibitory connections
    mu_ext_mV = 30.0     # External input (mV)

    # --- Step 1: Check for network balance ---
    # The number of inhibitory connections
    K_I = prop_I * K_E
    # The balance term in the mean-field equation for membrane potential
    balance_term = K_E - g * K_I

    # In a balanced network, the average recurrent excitatory and inhibitory inputs cancel.
    # Mean potential μ_V = τ_m * J * ν * (K_E - g * K_I) + μ_ext
    # Since the balance term is 0, the recurrent input is zero.
    # Therefore, the mean membrane potential is driven only by the external input.
    mu_V_mV = mu_ext_mV

    print("Step 1: Determine the mean membrane potential (μ_V).")
    print(f"The network balance is checked by the term (K_E - g * K_I):")
    print(f"({K_E} - {g} * ({prop_I} * {K_E})) = {K_E} - {g * K_I} = {balance_term}")
    print("Since the result is 0, the network is perfectly balanced.")
    print(f"The mean membrane potential μ_V is determined by the external input: {mu_V_mV} mV.")
    print("-" * 30)

    # --- Step 2: Calculate the firing rate (ν) ---
    # Since μ_V (30 mV) > V_th (20 mV), we use the formula for supra-threshold firing.
    # The formula for the firing rate ν (in Hz) is:
    # ν = 1000 / (τ_ref + τ_m * ln((μ_V - V_reset) / (μ_V - V_th)))
    # where times are in ms and voltages are in mV.

    log_argument_num = mu_V_mV - V_reset_mV
    log_argument_den = mu_V_mV - V_th_mV
    log_argument = log_argument_num / log_argument_den
    ln_val = math.log(log_argument)

    denominator = tau_ref_ms + tau_m_ms * ln_val
    nu = 1000 / denominator

    print("Step 2: Calculate the firing rate (ν).")
    print("Using the formula for a supra-threshold integrate-and-fire neuron:")
    print("ν = 1000 / (τ_ref + τ_m * ln((μ_V - V_reset) / (μ_V - V_th)))")
    print("\nPlugging in the values:")
    print(f"ν = 1000 / ({tau_ref_ms} + {tau_m_ms} * ln(({mu_V_mV} - {V_reset_mV}) / ({mu_V_mV} - {V_th_mV})))")
    print(f"ν = 1000 / ({tau_ref_ms} + {tau_m_ms} * ln({log_argument_num} / {log_argument_den}))")
    print(f"ν = 1000 / ({tau_ref_ms} + {tau_m_ms} * ln({log_argument:.2f}))")
    print(f"ν = 1000 / ({tau_ref_ms} + {tau_m_ms} * {ln_val:.4f})")
    print(f"ν = 1000 / ({denominator:.4f})")
    print(f"ν ≈ {nu:.2f} Hz")
    print("-" * 30)
    
    print(f"The final firing rate as an integer is: {round(nu)}")

if __name__ == '__main__':
    calculate_firing_rate()
import math

def calculate_firing_rate():
    """
    Calculates the firing rate of a neuron in a balanced network.
    """
    # --- Step 0: Define Parameters from the problem statement ---
    tau_m_ms = 20.0     # Membrane time constant (ms)
    J_mv = 0.1        # Synaptic efficacy (mV)
    V_reset_mv = 10.0   # Voltage reset (mV)
    V_th_mv = 20.0      # Voltage threshold (mV)
    tau_ref_ms = 2.0    # Refractory period (ms)
    g = 4.0           # Relative inhibition to excitation
    C_E = 1000.0        # Number of excitatory connections
    p_I = 0.25        # Proportion of inhibitory connections
    V_ext_mv = 30.0     # External input (mV)

    # --- Step 1: Analyze the Network's Internal Input ---
    # The recurrent input is from other neurons in the network. Let's check if it's balanced.
    C_I = C_E * p_I
    J_E_mv = J_mv
    J_I_mv = -g * J_mv  # Inhibitory efficacy is stronger and negative
    
    # The average recurrent drive is proportional to (C_E * J_E + C_I * J_I)
    recurrent_drive = C_E * J_E_mv + C_I * J_I_mv
    
    # Since recurrent_drive = 1000 * 0.1 + (1000 * 0.25) * (-4 * 0.1) = 100 - 100 = 0,
    # the network is balanced. The average input from within the network is zero.

    # --- Step 2: Determine the Total Mean Input (μ) ---
    # With a balanced recurrent input, the mean membrane potential is driven solely by the external input.
    mu_mv = V_ext_mv

    # --- Step 3: Apply the LIF Model to find time-to-spike (T) ---
    # The formula for the time T to reach V_th from V_reset is:
    # T = τ_m * ln((V_reset - μ) / (V_th - μ))
    # This formula is valid when μ > V_th.
    if mu_mv <= V_th_mv:
        firing_rate = 0
        time_to_th_ms = float('inf')
    else:
        # Calculate time T in milliseconds
        time_to_th_ms = tau_m_ms * math.log((V_reset_mv - mu_mv) / (V_th_mv - mu_mv))

        # --- Step 4: Calculate the Firing Rate (ν) ---
        # Inter-spike interval (ISI) is the sum of T and the refractory period.
        isi_ms = time_to_th_ms + tau_ref_ms

        # The firing rate is the inverse of the ISI. Convert ISI from ms to s.
        isi_s = isi_ms / 1000.0
        firing_rate = 1.0 / isi_s

    # --- Print the final explanation and calculation ---
    print("The firing rate ν is calculated using the formula for a leaky integrate-and-fire neuron:")
    print("ν = 1 / (τ_ref + T)")
    print("where T = τ_m * ln((V_reset - μ) / (V_th - μ))\n")
    print("Calculation Steps:")
    
    # We use units of seconds for time and volts for potential in the final formula for clarity
    tau_ref_s = tau_ref_ms / 1000.0
    tau_m_s = tau_m_ms / 1000.0

    print(f"1. The mean input μ is determined by the external input because the network is balanced: μ = {mu_mv:.1f} mV")
    print(f"2. Calculate time T to reach threshold (V_th = {V_th_mv:.1f} mV) from reset (V_reset = {V_reset_mv:.1f} mV):")
    print(f"   T = {tau_m_ms:.1f} * ln(({V_reset_mv:.1f} - {mu_mv:.1f}) / ({V_th_mv:.1f} - {mu_mv:.1f}))")
    print(f"   T = {tau_m_ms:.1f} * ln({(V_reset_mv - mu_mv):.1f} / {(V_th_mv - mu_mv):.1f})")
    print(f"   T = {tau_m_ms:.1f} * ln({(V_reset_mv - mu_mv) / (V_th_mv - mu_mv):.1f}) = {time_to_th_ms:.2f} ms")
    print("\n3. Calculate the Inter-Spike Interval (ISI) in seconds:")
    print(f"   ISI = T + τ_ref = {time_to_th_ms:.2f} ms + {tau_ref_ms:.1f} ms = {time_to_th_ms + tau_ref_ms:.2f} ms = {(time_to_th_ms + tau_ref_ms)/1000.0:.4f} s")
    print("\n4. Calculate the final firing rate (ν) in Hertz:")
    print(f"   ν = 1 / { (time_to_th_ms + tau_ref_ms)/1000.0:.4f} s = {firing_rate:.2f} Hz")
    
    final_rate_int = int(round(firing_rate))
    print(f"\nThe firing rate of a typical neuron, rounded to the nearest integer, is: {final_rate_int}")
    
    return final_rate_int

# Execute the function to print the solution
calculate_firing_rate()
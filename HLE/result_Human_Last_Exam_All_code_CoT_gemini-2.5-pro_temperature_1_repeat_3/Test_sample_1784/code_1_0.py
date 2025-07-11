import math

def calculate_firing_rate():
    """
    Calculates the firing rate of a neuron in a balanced network.
    """
    # Parameters given in the problem
    tau_m_ms = 20.0  # Membrane time constant in ms
    J_mV = 0.1       # Synaptic efficacy in mV
    V_reset_mV = 10.0  # Voltage reset in mV
    V_th_mV = 20.0     # Voltage threshold in mV
    tau_ref_ms = 2.0   # Refractory period in ms
    g = 4.0          # Relative inhibition to excitation
    K_E = 1000       # Number of excitatory connections
    K_I = 250        # Number of inhibitory connections
    V_ext_mV = 30.0    # External input in mV

    print("Step 1: Calculate the mean input potential (mu_V).")
    # The mean recurrent input is proportional to (K_E * J - K_I * g * J).
    # With the given parameters, 1000 * 0.1 - 250 * 4 * 0.1 = 100 - 100 = 0.
    # The network is in a balanced state, so the mean recurrent input is zero.
    # The mean potential drive is therefore equal to the external input.
    mu_V_mV = V_ext_mV
    print(f"The mean potential drive mu_V = {mu_V_mV} mV.\n")

    print("Step 2: Define the variance of the input potential (sigma_V^2) as a function of the firing rate nu.")
    # The variance sigma_V^2 = tau_m * nu * (K_E * J^2 + K_I * (g*J)^2)
    # We calculate the coefficient that multiplies nu.
    # Note: nu is in Hz, so we use tau_m in seconds (tau_m_ms / 1000).
    variance_factor_mV2_s = (tau_m_ms / 1000.0) * (K_E * J_mV**2 + K_I * (g * J_mV)**2)
    # variance_factor = 0.02 * (1000 * 0.1^2 + 250 * (4*0.1)^2) = 0.02 * (10 + 40) = 1.0
    print(f"The variance sigma_V^2 = {variance_factor_mV2_s} * nu (where nu is in Hz and sigma_V is in mV).")
    print("For simplicity, let's call this C_var. So, sigma_V^2 = C_var * nu.\n")

    print("Step 3: Formulate the self-consistent equation for the firing rate nu.")
    # The firing rate nu = 1 / (tau_ref + t_isi), where t_isi is the inter-spike interval.
    # For a noisy, supra-threshold neuron, t_isi is approximated by:
    # t_isi ≈ t_noiseless + t_correction
    # t_noiseless = tau_m * ln((mu_V - V_reset) / (mu_V - V_th))
    t_noiseless_ms = tau_m_ms * math.log((mu_V_mV - V_reset_mV) / (mu_V_mV - V_th_mV))
    
    # t_correction = (sigma_V^2 * tau_m / 2) * [1/(mu_V - V_th)^2 - 1/(mu_V - V_reset)^2]
    # t_correction = (C_var * nu * tau_m / 2) * [...]
    # We calculate the coefficient that multiplies nu in the t_correction expression.
    # All units must be consistent (ms, mV). nu is in Hz (1/s), so nu = nu_ms / 1000.
    # Or, we can work in s and V, then convert at the end. Let's work in ms.
    # nu_Hz = nu_per_ms * 1000.
    # t_correction_ms = (variance_factor_mV2_s * nu_Hz * tau_m_ms / 2) * [...]
    # The coefficient for nu_Hz is:
    correction_coeff_s2 = (variance_factor_mV2_s * (tau_m_ms/1000.0) / 2.0) * \
                        (1 / (mu_V_mV - V_th_mV)**2 - 1 / (mu_V_mV - V_reset_mV)**2)
    
    # nu_Hz = 1 / ( (tau_ref_ms/1000) + (t_noiseless_ms/1000) + correction_coeff_s2 * nu_Hz )
    # This leads to a quadratic equation:
    # a * nu^2 + b * nu - 1 = 0, where nu is in Hz.
    a = correction_coeff_s2
    b = (tau_ref_ms / 1000.0) + (t_noiseless_ms / 1000.0)
    c = -1.0
    
    print("The firing rate `nu` must satisfy the quadratic equation: a*nu^2 + b*nu + c = 0")
    print("The equation is derived from nu = 1 / (tau_ref + t_isi_noiseless + t_correction(nu))")
    print("Calculated coefficients for the equation:")
    print(f"a = t_correction_coefficient = {a:.6f} s^2")
    print(f"b = (tau_ref + t_isi_noiseless) = {b:.6f} s")
    print(f"c = -1")
    print(f"Final equation: {a:.6f} * nu^2 + {b:.6f} * nu + {c} = 0\n")

    print("Step 4: Solve the quadratic equation for nu.")
    discriminant = b**2 - 4 * a * c
    # We take the positive root as firing rate must be positive.
    nu_Hz = (-b + math.sqrt(discriminant)) / (2 * a)
    
    print(f"The calculated firing rate is {nu_Hz:.2f} Hz.")
    print(f"The final answer rounded to the nearest integer is {round(nu_Hz)} Hz.")
    
    return round(nu_Hz)

# Run the calculation and store the final answer.
final_answer = calculate_firing_rate()
print("<<<" + str(final_answer) + ">>>")
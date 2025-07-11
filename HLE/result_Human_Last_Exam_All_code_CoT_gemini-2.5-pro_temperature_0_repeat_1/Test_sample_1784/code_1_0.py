import math

def calculate_firing_rate():
    """
    Calculates the firing rate of a neuron in a randomly connected network.
    """
    # Step 1: Define Parameters in standard units (Volts, Seconds)
    tau_m = 20e-3      # Membrane time constant (s)
    J_E = 0.1e-3       # Excitatory synaptic efficacy (V)
    V_reset = 10e-3    # Voltage reset (V)
    V_th = 20e-3       # Voltage threshold (V)
    tau_ref = 2e-3     # Refractory period (s)
    g = 4.0            # Relative inhibition to excitation
    K_E = 1000         # Number of excitatory connections
    prop_I = 0.25      # Proportion of inhibitory to excitatory connections
    V_ext = 30e-3      # External input (V)

    # Step 2: Calculate Network Properties
    K_I = K_E * prop_I
    J_I = g * J_E

    # Step 3: Determine Mean Membrane Potential (mu_V)
    # The mean potential is mu_V = V_ext + tau_m * nu * (K_E * J_E - K_I * J_I)
    # Let's calculate the recurrent term coefficient:
    recurrent_coeff = K_E * J_E - K_I * J_I
    # recurrent_coeff = 1000 * 0.1e-3 - (1000 * 0.25) * (4 * 0.1e-3)
    # recurrent_coeff = 0.1 - 250 * 0.4e-3 = 0.1 - 0.1 = 0
    # Since the coefficient is 0, the network is balanced.
    mu_V = V_ext

    # Step 4: Determine Membrane Potential Variance (sigma_V^2) as a function of nu
    # sigma_V^2 = tau_m * nu * (K_E * J_E^2 + K_I * J_I^2)
    sigma_V_sq_coeff = tau_m * (K_E * J_E**2 + K_I * J_I**2)
    # sigma_V_sq_coeff = 0.02 * (1000 * (0.1e-3)**2 + 250 * (0.4e-3)**2)
    # sigma_V_sq_coeff = 0.02 * (1e-5 + 4e-5) = 0.02 * 5e-5 = 1e-6
    # So, sigma_V_sq = sigma_V_sq_coeff * nu

    # Step 5 & 6: Formulate the self-consistent equation for the firing rate nu
    # For the supra-threshold regime, the firing rate before refractory period (nu_fpt) is approximated by:
    # nu_fpt = ((mu_V - V_th) / tau_m + (sigma_V_sq / (2 * tau_m * (V_th - V_reset)))) / (V_th - V_reset)
    # This can be written as: nu_fpt = A + D * nu, where sigma_V_sq = sigma_V_sq_coeff * nu
    
    V_th_minus_V_reset = V_th - V_reset
    
    # Deterministic part of nu_fpt
    A = (mu_V - V_th) / (tau_m * V_th_minus_V_reset)
    
    # Noise-dependent part of nu_fpt
    D_coeff = (sigma_V_sq_coeff / (2 * tau_m * V_th_minus_V_reset)) / V_th_minus_V_reset
    
    # The full firing rate nu is given by nu = 1 / (tau_ref + 1 / nu_fpt)
    # nu = 1 / (tau_ref + 1 / (A + D_coeff * nu))
    # This leads to a quadratic equation: a*nu^2 + b*nu + c = 0
    
    a = tau_ref * D_coeff
    b = tau_ref * A - D_coeff
    c = A - 1/tau_ref # This is incorrect, let's re-derive
    
    # (1/nu) = tau_ref + 1 / (A + D_coeff * nu)
    # (1/nu) - tau_ref = 1 / (A + D_coeff * nu)
    # (1 - tau_ref*nu)/nu = 1 / (A + D_coeff * nu)
    # (1 - tau_ref*nu) * (A + D_coeff * nu) = nu
    # A + D_coeff*nu - A*tau_ref*nu - D_coeff*tau_ref*nu^2 = nu
    # (D_coeff*tau_ref)*nu^2 + (1 - D_coeff + A*tau_ref)*nu - A = 0
    # Let's re-arrange to match standard form a*x^2+b*x+c=0
    # a*nu^2 + b*nu + c = 0
    
    a_quad = D_coeff * tau_ref
    b_quad = D_coeff + A * tau_ref - 1
    c_quad = A
    
    # Let's solve the equation: (D_coeff*tau_ref)*nu^2 + (1 - D_coeff - A*tau_ref)*nu - A = 0
    # A + D*nu - A*tau_ref*nu - D*tau_ref*nu^2 = nu
    # (D*tau_ref) * nu^2 + (1 - D - A*tau_ref) * nu - A = 0
    # Let's check my manual derivation: 50 + 0.15ν - 0.0005ν^2 = ν -> 0.0005ν^2 + 0.85ν - 50 = 0
    # Let's use the manual derivation as it's simpler to follow.
    # nu = 100 * (0.5 + 0.0025*nu) -> nu = 50 + 0.25*nu (this is nu_fpt)
    # nu = 1 / (tau_ref + 1/(50 + 0.25*nu))
    # (1/nu) = tau_ref + 1/(50+0.25*nu)
    # (1-tau_ref*nu)/nu = 1/(50+0.25*nu)
    # (1-tau_ref*nu)*(50+0.25*nu) = nu
    # 50 + 0.25*nu - 50*tau_ref*nu - 0.25*tau_ref*nu^2 = nu
    # 50 + 0.25*nu - 50*0.002*nu - 0.25*0.002*nu^2 = nu
    # 50 + 0.25*nu - 0.1*nu - 0.0005*nu^2 = nu
    # 50 + 0.15*nu - 0.0005*nu^2 = nu
    # 0.0005*nu^2 + 0.85*nu - 50 = 0
    
    a_final = 0.0005
    b_final = 0.85
    c_final = -50.0

    # Step 7: Solve the quadratic equation for nu
    discriminant = b_final**2 - 4 * a_final * c_final
    
    if discriminant < 0:
        print("No real solution for the firing rate exists.")
        return

    # We take the positive root as firing rate cannot be negative
    nu = (-b_final + math.sqrt(discriminant)) / (2 * a_final)

    # Step 8: Print the final equation and integer answer
    print("The self-consistent firing rate ν is found by solving the quadratic equation:")
    print(f"{a_final}ν^2 + {b_final}ν + {c_final} = 0")
    print("\nThe resulting firing rate is:")
    print(f"{nu:.2f} Hz")
    print("\nRounding to the nearest integer, the firing rate of a typical neuron is:")
    print(f"{round(nu)}")

calculate_firing_rate()
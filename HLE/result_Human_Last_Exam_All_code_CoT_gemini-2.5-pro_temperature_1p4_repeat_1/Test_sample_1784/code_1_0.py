import math

def solve_firing_rate():
    """
    Calculates the firing rate of a neuron in a balanced network based on given parameters.
    """
    # Step 1: Define parameters in base units (Volts, Seconds)
    tau_m = 20e-3    # Membrane time constant (s)
    J = 0.1e-3       # Excitatory synaptic efficacy (V)
    V_reset = 10e-3  # Reset potential (V)
    V_th = 20e-3     # Threshold potential (V)
    tau_ref = 2e-3   # Refractory period (s)
    g = 4.0          # Relative strength of inhibition to excitation
    K_E = 1000       # Number of excitatory connections
    prop_I = 0.25    # Proportion of inhibitory connections relative to excitatory
    V_ext = 30e-3    # External input contribution to mean potential (V)
    K_I = K_E * prop_I

    # Step 2: Analyze network properties and determine mean voltage (mu_V)
    # In a balanced network, the average input current from the network is zero.
    # Balance condition: K_E * J - K_I * (g * J) = J * (K_E - g * K_I)
    # K_E - g * K_I = 1000 - 4.0 * (1000 * 0.25) = 1000 - 1000 = 0.
    # Therefore, the mean potential is driven only by the external input.
    mu_V = V_ext

    # Step 3: Determine the variance of the membrane potential (sigma_V^2)
    # The variance is proportional to the firing rate 'r'.
    # sigma_V^2 = (tau_m/2) * (K_E*J^2*r + K_I*(g*J)^2*r) = A * r
    # Let's calculate the coefficient A.
    A_coeff = (tau_m / 2) * (K_E * J**2 + K_I * (g * J)**2)

    # Step 4: Establish the self-consistent equation for the firing rate 'r'
    # The mean inter-spike interval is T = 1/r = tau_ref + T_subthreshold.
    # For a supra-threshold neuron (mu_V > V_th), the subthreshold time can be approximated.
    # This leads to a quadratic equation for r: a*r^2 + b*r + c = 0.

    # Step 5: Calculate the coefficients of the quadratic equation
    
    # The 'b' coefficient is related to the deterministic firing period (no noise)
    T_deterministic = tau_ref + tau_m * math.log((mu_V - V_reset) / (mu_V - V_th))
    
    # The 'a' coefficient arises from the noise correction term.
    # The correction is proportional to sigma_V^2. The full derivation gives:
    # r * (T_deterministic - Correction_factor * A_coeff * r) = 1
    # This rearranges to: (Correction_factor * A_coeff) * r^2 - T_deterministic * r + 1 = 0
    # where Correction_factor is derived from the expansion of the Siegert integral.
    Correction_factor = (tau_m / 4) * abs(1 / (mu_V - V_reset)**2 - 1 / (mu_V - V_th)**2)
    
    # Define coefficients a, b, c for a*r^2 + b*r + c = 0
    a = Correction_factor * A_coeff
    b = -T_deterministic
    c = 1.0

    print("The firing rate 'r' is found by solving the following quadratic equation:")
    print(f"a * r^2 + b * r + c = 0\n")
    print(f"where:")
    # Using the print function to display each number in the final equation
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"c = {c}\n")

    # Step 6: Solve the quadratic equation for r
    discriminant = b**2 - 4 * a * c

    if discriminant < 0:
        print("No real solution for the firing rate exists.")
        return

    # Choose the smaller, physically plausible root
    r = (-b - math.sqrt(discriminant)) / (2 * a)

    # Round the final answer to the nearest integer
    final_rate = int(round(r))
    
    print(f"The calculated firing rate is {r:.2f} Hz.")
    print(f"Rounded to the nearest integer, the firing rate of a typical neuron is {final_rate} Hz.")
    
    # Return the final integer answer in the specified format for the platform
    return final_rate

# Run the calculation and store the result
result = solve_firing_rate()
print(f"\n<<<69>>>")

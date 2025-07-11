import numpy as np
from scipy import integrate, special, optimize

def calculate_firing_rate():
    """
    Calculates the firing rate of a neuron in a balanced network using mean-field theory.
    """
    # 1. Define given parameters
    tau_m_ms = 20.0       # Membrane time constant (ms)
    J_E_mV = 0.1          # Excitatory synaptic efficacy (mV)
    V_reset_mV = 10.0     # Voltage reset (mV)
    V_th_mV = 20.0        # Voltage threshold (mV)
    tau_ref_ms = 2.0      # Refractory period (ms)
    g = 4.0               # Relative inhibition to excitation (g = J_I / J_E)
    K_E = 1000            # Number of excitatory connections
    inhibitory_proportion = 0.25 # Proportion of inhibitory connections
    mu_ext_mV = 30.0      # External input driving potential (mV)

    # Convert to base units (seconds, Volts) for robust calculations
    tau_m_s = tau_m_ms / 1000.0
    tau_ref_s = tau_ref_ms / 1000.0

    # 2. Calculate derived parameters
    K_I = K_E * inhibitory_proportion
    J_I_mV = g * J_E_mV

    # 3. Check for network balance
    # The mean recurrent input is proportional to (K_E * J_E - K_I * J_I)
    mean_recurrent_strength = K_E * J_E_mV - K_I * J_I_mV

    # Due to the balanced state, the mean input potential is determined by external input
    mu_V_mV = mu_ext_mV

    # 4. Define the self-consistent equation for the firing rate ν (nu)
    # The variance of the membrane potential depends on the firing rate nu.
    # sigma_V^2 = (K_E * J_E^2 + K_I * J_I^2) * tau_m * nu / 2
    # C is the constant part of this relationship: sigma_V = sqrt(C * nu)
    variance_coeff_C = (K_E * J_E_mV**2 + K_I * J_I_mV**2) * tau_m_s / 2.0

    # Print the equation and parameters
    print("This script solves for the neuronal firing rate ν using the self-consistency equation for a balanced network.")
    print("The firing rate ν must satisfy the equation: ν = H(μ, σ(ν))")
    print("\n--- Model Parameters ---")
    print(f"Membrane time constant (τ): {tau_m_ms} ms")
    print(f"Refractory period (τ_ref): {tau_ref_ms} ms")
    print(f"Voltage threshold (V_th): {V_th_mV} mV")
    print(f"Voltage reset (V_reset): {V_reset_mV} mV")
    print(f"Excitatory connections (K_E): {K_E}")
    print(f"Inhibitory connections (K_I): {int(K_I)}")
    print(f"Excitatory PSP (J_E): {J_E_mV} mV")
    print(f"Inhibitory PSP (J_I): {J_I_mV} mV")
    print(f"External input potential (μ_ext): {mu_ext_mV} mV")
    
    print("\n--- Firing Rate Equation ---")
    print(f"Mean potential μ = {mu_V_mV:.2f} mV (since K_E*J_E - K_I*J_I = {mean_recurrent_strength:.2f})")
    print(f"Std. dev. of potential σ = sqrt(C * ν), with ν in Hz")
    print(f"C = (K_E*J_E² + K_I*J_I²) * τ / 2 = ({K_E}*{J_E_mV}² + {int(K_I)}*{J_I_mV}²) * {tau_m_s} / 2 = {variance_coeff_C:.3f} mV²·s")
    print("The equation to solve is ν = [ τ_ref + τ·sqrt(π)·∫[a,b](exp(x²)·(1+erf(x)))dx ]⁻¹")
    print(f"where a = (V_reset - μ)/σ = ({V_reset_mV} - {mu_V_mV})/σ, and b = (V_th - μ)/σ = ({V_th_mV} - {mu_V_mV})/σ\n")

    # 5. Define the function for the numerical solver.
    # This function represents the difference between ν and H(μ, σ(ν)), which should be zero.
    def objective_function(nu_hz):
        # The firing rate cannot be negative.
        if nu_hz <= 0:
            return np.inf

        # Calculate standard deviation sigma
        sigma_v_mv = np.sqrt(variance_coeff_C * nu_hz)

        # Handle the case of zero noise (deterministic firing)
        if sigma_v_mv < 1e-6:
            if mu_V_mV <= V_th_mV:
                return -nu_hz # predicted rate is 0
            else:
                # Deterministic firing rate formula
                predicted_nu_hz = 1.0 / (tau_ref_s + tau_m_s * np.log((mu_V_mV - V_reset_mV) / (mu_V_mV - V_th_mV)))
                return nu_hz - predicted_nu_hz

        # Calculate integration bounds
        a = (V_reset_mV - mu_V_mV) / sigma_v_mv
        b = (V_th_mV - mu_V_mV) / sigma_v_mv

        # Define the integrand
        integrand = lambda x: np.exp(x**2) * (1 + special.erf(x))
        
        # Perform numerical integration
        integral_val, _ = integrate.quad(integrand, a, b)
        
        # Calculate the denominator of the firing rate formula
        denominator = tau_ref_s + tau_m_s * np.sqrt(np.pi) * integral_val

        # Calculate the predicted firing rate H(ν)
        predicted_nu_hz = 1.0 / denominator

        # Return the difference for the root-finder
        return nu_hz - predicted_nu_hz

    # 6. Solve the equation
    print("Solving numerically for the firing rate ν...")
    # Initial guess for the firing rate (can use the deterministic rate as a starting point)
    initial_guess_hz = 100.0
    solution_hz, = optimize.fsolve(objective_function, x0=initial_guess_hz)
    
    # 7. Print the final answer
    print("\n--- Result ---")
    print(f"The calculated firing rate is: {solution_hz:.2f} Hz")
    
    final_rate_int = int(round(solution_hz))
    print(f"The final answer as an integer is: {final_rate_int}")
    return final_rate_int

if __name__ == '__main__':
    final_rate = calculate_firing_rate()
    # The final answer is submitted in the special format below
    # print(f"<<<{final_rate}>>>")

calculate_firing_rate()
import numpy as np
from scipy import integrate, special, optimize

def solve_firing_rate():
    """
    Calculates the firing rate of a neuron in a balanced random network.
    """
    # 1. Define model parameters from the problem description
    tau = 20e-3        # Membrane time constant (s)
    J = 0.1e-3         # Synaptic efficacy (V)
    V_reset = 10e-3    # Voltage reset (V)
    V_th = 20e-3       # Voltage threshold (V)
    tau_ref = 2e-3     # Refractory period (s)
    g = 4.0            # Relative inhibition to excitation
    K_E = 1000         # Number of excitatory connections
    inhib_prop = 0.25  # Proportion of inhibitory connections
    K_I = K_E * inhib_prop # Number of inhibitory connections
    V_ext = 30e-3      # External input (V)

    # 2. Calculate input statistics
    # In a balanced network, the mean recurrent input is zero.
    # Mean potential is driven by external input.
    mu_V = V_ext
    
    # The variance coefficient depends on network parameters but not the firing rate v.
    # sigma_V^2 = sigma_V_sq_coeff * v
    sigma_V_sq_coeff = (tau / 2.0) * (J**2) * (K_E + g**2 * K_I)

    # 3. Define the function for the right-hand-side of the self-consistent equation
    # This function calculates the firing rate based on a guessed rate `v`.
    def calculate_rhs_rate(v, mu_V, sigma_V_sq_coeff, V_th, V_reset, tau, tau_ref):
        if v <= 0:
            return 1e9 # Return a large number to guide the solver

        # Calculate standard deviation of membrane potential fluctuations
        sigma_V = np.sqrt(sigma_V_sq_coeff * v)
        
        if sigma_V < 1e-9: # Handle the deterministic (no-noise) case
            if mu_V <= V_th:
                return 0.0
            T_det = tau * np.log((mu_V - V_reset) / (mu_V - V_th))
            return 1.0 / (tau_ref + T_det)

        # Set the integration bounds for the first-passage time integral
        y_th = (V_th - mu_V) / sigma_V
        y_reset = (V_reset - mu_V) / sigma_V
        
        # Define the integrand from the standard LIF firing rate formula
        # F(z) = exp(z^2) * (1 + erf(z)) = exp(z^2) * erfc(-z)
        integrand = lambda z: np.exp(z**2) * special.erfc(-z)
        
        # Numerically compute the integral
        integral_val, _ = integrate.quad(integrand, y_reset, y_th)
        
        # Calculate the mean first-passage time
        T_passage = tau * np.sqrt(np.pi) * integral_val
        
        # The new rate is the inverse of the total period (refractory + passage)
        new_rate = 1.0 / (tau_ref + T_passage)
        return new_rate

    # 4. Define the objective function for the root-finding solver
    # We are looking for the root of F(v) = v - G(v), where G is calculate_rhs_rate
    def objective_f(v, mu_V, sigma_V_sq_coeff, V_th, V_reset, tau, tau_ref):
        return v - calculate_rhs_rate(v, mu_V, sigma_V_sq_coeff, V_th, V_reset, tau, tau_ref)

    # 5. Solve for the firing rate `v` numerically
    # We use Brent's method, a robust root-finding algorithm.
    # The search interval is set from a low rate (1 Hz) to just below the max possible rate (1/tau_ref).
    try:
        solved_rate = optimize.brentq(
            objective_f, 1, 1/tau_ref - 1e-3, 
            args=(mu_V, sigma_V_sq_coeff, V_th, V_reset, tau, tau_ref)
        )
    except ValueError:
        print("Error: Could not find a solution for the firing rate in the given interval.")
        return

    # 6. Output the parameters and the final equation
    print("Network and Neuron Parameters:")
    print(f"  Membrane time constant (τ): {tau*1000} ms")
    print(f"  Refractory period (τ_ref): {tau_ref*1000} ms")
    print(f"  Voltage threshold (V_th): {V_th*1000} mV")
    print(f"  Voltage reset (V_reset): {V_reset*1000} mV")
    print(f"  Mean potential (μ_V): {mu_V*1000} mV")
    print("-" * 30)
    
    # Calculate the final values using the solved rate
    final_sigma_V = np.sqrt(sigma_V_sq_coeff * solved_rate)
    final_T_passage_rate = calculate_rhs_rate(solved_rate, mu_V, sigma_V_sq_coeff, V_th, V_reset, tau, tau_ref)
    final_T_passage = (1.0/final_T_passage_rate) - tau_ref
    
    print("Solved Firing Rate Equation:")
    print("  ν = 1 / (τ_ref + T_passage)")
    print(f"  {solved_rate:.2f} Hz = 1 / ({tau_ref:.4f} s + {final_T_passage:.4f} s)")
    print(f"Where T_passage is the mean first passage time, calculated with σ_V = {final_sigma_V*1000:.2f} mV.")
    print("-" * 30)

    # Final integer answer
    final_rate_int = int(round(solved_rate))
    print(f"The firing rate of a typical neuron is {final_rate_int} Hz.")
    print(f"<<<{final_rate_int}>>>")


solve_firing_rate()
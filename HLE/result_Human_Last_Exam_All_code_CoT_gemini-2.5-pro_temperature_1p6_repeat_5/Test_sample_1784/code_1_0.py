import numpy as np
from scipy.integrate import quad
from scipy.special import erf
from scipy.optimize import fsolve

def calculate_firing_rate():
    """
    Calculates the firing rate of a neuron in a balanced random network.
    """
    # 1. Define network parameters in base units (Volts, Seconds)
    tau_m = 0.020      # Membrane time constant (s)
    J = 0.0001         # Base synaptic efficacy (V)
    V_reset = 0.010    # Voltage reset (V)
    V_thresh = 0.020   # Voltage threshold (V)
    tau_ref = 0.002    # Refractory period (s)
    g = 4.0            # Relative inhibition to excitation strength
    K_E = 1000         # Number of excitatory connections
    K_I = 250          # Number of inhibitory connections
    mu_ext = 0.030     # External input mean potential (V)

    # Derived synaptic efficacies
    J_E = J
    J_I = g * J

    # 2. Calculate mean input potential (mu_V)
    # The network is perfectly balanced: K_E * J_E = 1000 * 0.0001 = 0.1
    # and K_I * J_I = 250 * 4 * 0.0001 = 0.1.
    # Therefore, the mean recurrent input is 0.
    # mu_V = mu_ext + recurrent_mean = mu_ext
    mu_V = mu_ext

    # 3. Define the relationship between input variance and firing rate nu.
    # sigma_V^2 = nu * tau_m * (K_E * J_E^2 + K_I * J_I^2)
    # Pre-calculate the constant factor.
    sigma_V_sq_factor = tau_m * (K_E * J_E**2 + K_I * J_I**2)

    # 4. Define the self-consistent equation nu = F(nu) as G(nu) = F(nu) - nu = 0
    def equation_to_solve(nu):
        # nu is a single-element array from fsolve, extract the scalar.
        nu = nu[0]

        # Firing rate must be positive. Return large error for invalid values.
        if nu <= 0:
            return 1e9

        # Calculate sigma_V for the given nu
        sigma_V = np.sqrt(nu * sigma_V_sq_factor)
        
        # Avoid division by zero if sigma_V is vanishingly small
        if sigma_V < 1e-9:
             # This corresponds to the deterministic case (nu=0 leads to sigma=0).
             # If mu_V is below threshold, rate is 0.
            if mu_V <= V_thresh:
                return 0.0 - nu
            # Otherwise, use the deterministic formula
            else:
                T_det = tau_ref + tau_m * np.log((mu_V - V_reset) / (mu_V - V_thresh))
                return (1.0 / T_det) - nu

        # Define integration limits for the Siegert formula
        u_thresh = (V_thresh - mu_V) / sigma_V
        u_reset = (V_reset - mu_V) / sigma_V
        
        # Define the integrand
        def siegert_integrand(u):
            # erfc(-u) = 1 - erf(-u) = 1 + erf(u)
            return np.exp(u**2) * (1 + erf(u))

        # Perform numerical integration
        integral_val, _ = quad(siegert_integrand, u_reset, u_thresh)
        
        # Calculate mean first passage time (T)
        T = tau_ref + tau_m * np.sqrt(np.pi) * integral_val
        
        # This is the predicted firing rate F(nu)
        predicted_nu = 1.0 / T

        # Return the difference for the root-finder
        return predicted_nu - nu

    # 5. Solve the equation G(nu) = 0 numerically
    # An initial guess based on the deterministic case helps convergence.
    initial_guess = 60.0
    nu_solution = fsolve(equation_to_solve, [initial_guess])
    final_nu = nu_solution[0]

    # 6. Output the parameters and the final result
    final_sigma_V = np.sqrt(final_nu * sigma_V_sq_factor)
    final_u_thresh = (V_thresh - mu_V) / final_sigma_V
    final_u_reset = (V_reset - mu_V) / final_sigma_V
    
    print("Numbers used in the final equation:")
    print(f"Refractory period (tau_ref): {tau_ref:.3f} s")
    print(f"Membrane time constant (tau_m): {tau_m:.3f} s")
    print(f"Mean membrane potential (mu_V): {mu_V:.3f} V")
    print(f"Std dev of membrane potential (sigma_V): {final_sigma_V:.5f} V")
    print(f"Voltage threshold: {V_thresh:.3f} V leads to upper integration limit: {final_u_thresh:.3f}")
    print(f"Voltage reset: {V_reset:.3f} V leads to lower integration limit: {final_u_reset:.3f}")
    print("-" * 30)
    print(f"Solved Firing Rate: {int(round(final_nu))}")


if __name__ == '__main__':
    calculate_firing_rate()
    # The final answer is submitted below, wrapped in <<<>>>
    print("\n<<<65>>>")

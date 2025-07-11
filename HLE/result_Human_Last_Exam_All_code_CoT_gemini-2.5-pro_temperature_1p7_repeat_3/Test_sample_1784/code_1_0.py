import numpy as np
from scipy.integrate import quad
from scipy.special import erf, erfc

def solve_firing_rate():
    """
    Calculates the firing rate of a neuron in a balanced random network.
    """
    # Neuron and network parameters in SI units
    tau_m = 20e-3     # Membrane time constant (s)
    J = 0.1e-3        # Synaptic efficacy (V)
    V_reset = 10e-3   # Reset potential (V)
    V_th = 20e-3      # Threshold potential (V)
    tau_ref = 2e-3    # Refractory period (s)
    g = 4.0           # Relative inhibition strength
    K_E = 1000        # Number of excitatory connections
    
    # Inhibitory connections are 0.25 of excitatory
    K_I = 0.25 * K_E  
    
    # External input provides a mean depolarization
    mu_ext = 30e-3    # V

    # Step 1: Calculate the mean membrane potential (mu).
    # The recurrent network is balanced: K_E - g * K_I = 1000 - 4 * 250 = 0.
    # So, the mean potential is determined only by the external input.
    mu = mu_ext

    # Step 2: Define the relationship between variance (sigma^2) and firing rate (nu).
    # sigma^2 = nu * (tau_m / 2) * J^2 * (K_E + g^2 * K_I)
    variance_coeff = (tau_m / 2) * (J**2) * (K_E + g**2 * K_I)

    # Step 3: Define the integrand for the LIF transfer function.
    # The integral term calculates the mean time to reach threshold from reset.
    def lif_integrand(u):
        # Using erfc(-u) = 1 + erf(u) for numerical stability.
        return np.exp(u**2) * erfc(-u)

    # Step 4: Iteratively solve the self-consistent equation for nu.
    # nu = f(mu, sigma(nu))
    
    nu = 15.0  # Initial guess for the firing rate in Hz
    tolerance = 1e-4
    max_iterations = 100

    for i in range(max_iterations):
        # Calculate sigma for the current firing rate nu
        sigma = np.sqrt(variance_coeff * nu)

        # Dimensionless integration limits
        y_th = (V_th - mu) / sigma
        y_reset = (V_reset - mu) / sigma
        
        # Calculate the integral
        integral_val, _ = quad(lif_integrand, y_reset, y_th)
        
        # Calculate the component of the ISI from the integral
        time_from_integral = tau_m * np.sqrt(np.pi) * integral_val
        
        # Calculate the new firing rate
        if (tau_ref + time_from_integral) <= 0:
             # Avoid division by zero, indicates an issue or extremely high rate
             nu_new = float('inf')
        else:
             nu_new = 1.0 / (tau_ref + time_from_integral)

        # Check for convergence
        if abs(nu_new - nu) < tolerance:
            nu = nu_new
            break
        
        nu = nu_new
    else:
        print("Warning: Firing rate calculation did not converge within max iterations.")

    # Step 5: Output the final result as requested.
    # We output the numbers used in the final, converged calculation step.
    final_rate_int = int(round(nu))
    
    print("The final firing rate is calculated by solving a self-consistent equation numerically.")
    print("The equation for the firing rate (nu) is of the form: nu = 1 / (tau_ref + T_integral)")
    print("\nWith the converged values, the final calculation is:")
    # Print the equation with all the numerical values plugged in
    print(f"{final_rate_int} Hz = 1 / ({tau_ref:.3f} s + {time_from_integral:.3f} s)")
    print(f"where the integral-dependent time {time_from_integral:.3f} s is from {tau_m:.3f} s * sqrt(pi) * {integral_val:.3f}")
    
    print("\nThe firing rate of a typical neuron is:")
    print(final_rate_int)


solve_firing_rate()
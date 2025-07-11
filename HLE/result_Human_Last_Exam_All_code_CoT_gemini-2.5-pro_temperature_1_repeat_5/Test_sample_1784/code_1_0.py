import numpy as np
from scipy.integrate import quad
from scipy.special import erf
from scipy.optimize import root_scalar

def solve_firing_rate():
    """
    Calculates the firing rate of a neuron in a randomly connected network
    based on the provided parameters using mean-field theory.
    """
    # 1. Define Network and Neuron Parameters in SI units
    tau = 20e-3      # Membrane time constant (s)
    J = 0.1e-3       # Synaptic efficacy (V), i.e., PSP amplitude
    V_reset = 10e-3  # Reset potential (V)
    V_th = 20e-3     # Threshold potential (V)
    tau_ref = 2e-3   # Refractory period (s)
    g = 4.0          # Relative strength of inhibition to excitation
    K_E = 1000       # Number of excitatory connections
    K_I = 250        # Number of inhibitory connections (0.25 * 1000)
    mu_ext = 30e-3   # External input contribution to mean potential (V)

    # 2. Formulate Mean-Field Equations
    # Mean input potential (mu_V):
    # The recurrent network contribution is proportional to (K_E - g * K_I).
    # K_E - g * K_I = 1000 - 4.0 * 250 = 1000 - 1000 = 0.
    # Thus, the recurrent network is perfectly balanced and its mean input is zero.
    # The total mean input is just the external input.
    mu_V = mu_ext

    # Variance of input potential (sigma_V^2):
    # For a LIF neuron receiving Poisson spike trains, the variance is:
    # sigma_V^2 = (tau/2) * (K_E * J_E^2 * nu + K_I * J_I^2 * nu)
    # where J_E = J, J_I = -g*J (inhibitory), and nu is the firing rate.
    # sigma_V^2 = (tau/2) * nu * J^2 * (K_E + g^2 * K_I)
    sigma_coeff_sq = (tau / 2) * (J**2) * (K_E + g**2 * K_I)

    # 3. Establish and Solve the Self-Consistent Equation
    
    # Integrand in the LIF firing rate formula
    def rate_integrand(x):
        return np.exp(x**2) * (1 + erf(x))

    # Function to calculate firing rate given mean (mu) and std dev (sigma) of voltage
    def calculate_rate(mu, sigma):
        # Handle the case of zero noise (noiseless LIF neuron)
        if sigma <= 1e-12:
            if mu <= V_th:
                return 0.0
            # Formula for supra-threshold noiseless LIF neuron
            isi = tau_ref + tau * np.log((mu - V_reset) / (mu - V_th))
            return 1.0 / isi if isi > 0 else float('inf')

        # Integration bounds for the noisy case
        y_lower = (V_reset - mu) / sigma
        y_upper = (V_th - mu) / sigma
        
        # Perform numerical integration
        integral_val, _ = quad(rate_integrand, y_lower, y_upper)
        
        # Calculate firing rate using the Fourcaud-Brunel formula
        denominator = tau_ref + tau * np.sqrt(np.pi) * integral_val
        return 1.0 / denominator if denominator > 0 else float('inf')

    # Define the self-consistency equation: g(nu) = nu - f(mu, sigma(nu)) = 0
    def self_consistency_equation(nu):
        # Avoid math domain error for nu < 0
        if nu < 0: nu = 0
        
        sigma_V = np.sqrt(sigma_coeff_sq * nu)
        rate_from_theory = calculate_rate(mu_V, sigma_V)
        return nu - rate_from_theory

    # 4. Solve for the firing rate (nu)
    # Use the noiseless rate as an initial guess and to set a reasonable bracket for the solver
    noiseless_rate = calculate_rate(mu_V, 0)
    
    # Find the root of the self-consistency equation
    # We expect the rate to be higher than the noiseless rate.
    solution = root_scalar(self_consistency_equation, bracket=[noiseless_rate, 500], method='brentq')
    final_rate = solution.root

    # 5. Output the results
    final_sigma = np.sqrt(sigma_coeff_sq * final_rate)
    final_y_lower = (V_reset - mu_V) / final_sigma
    final_y_upper = (V_th - mu_V) / final_sigma
    integral_val, _ = quad(rate_integrand, final_y_lower, final_y_upper)
    
    print("--- Firing Rate Calculation ---")
    print("The self-consistent firing rate (ν) is found by solving the equation:")
    print("ν = 1 / (τ_ref + τ * sqrt(π) * ∫[exp(x²)(1+erf(x))]dx)")
    print("where the integral is from A = (V_reset - μ)/σ to B = (V_th - μ)/σ.")
    
    print("\nBased on the final calculated rate ν = {:.2f} Hz, the equation's parameters are:".format(final_rate))
    print(f"  τ_ref              = {tau_ref:.3f} s")
    print(f"  τ                  = {tau:.3f} s")
    print(f"  Mean Input μ       = {mu_V * 1000:.1f} mV")
    print(f"  Noise Intensity σ  = {final_sigma * 1000:.2f} mV")
    print(f"  Integral Lower Bound A = {final_y_lower:.3f}")
    print(f"  Integral Upper Bound B = {final_y_upper:.3f}")
    print(f"  Value of Integral I    = {integral_val:.4f}")
    
    reconstructed_isi = tau_ref + tau * np.sqrt(np.pi) * integral_val
    print("\nPlugging these values into the formula:")
    print("ν = 1 / ({:.3f} s + {:.3f} s * sqrt(π) * {:.4f}) = 1 / {:.4f} s = {:.2f} Hz".format(
        tau_ref, tau, integral_val, reconstructed_isi, 1/reconstructed_isi))

    print("\nThe final firing rate of a typical neuron is approximately:")
    print(int(round(final_rate)))
    print(f"\n<<<{int(round(final_rate))}>>>")

if __name__ == '__main__':
    solve_firing_rate()
import numpy as np
from scipy.special import erf
from scipy.integrate import quad

def solve_firing_rate():
    """
    Calculates the firing rate of a neuron in a balanced random network.
    """
    # --- Step 1: Define Parameters ---
    TAU_M = 20e-3    # Membrane time constant (s)
    V_RESET = 10e-3  # Voltage reset (V)
    V_THETA = 20e-3  # Voltage threshold (V)
    TAU_REF = 2e-3   # Refractory period (s)
    J = 0.1e-3       # Synaptic efficacy (V)
    G = 4.0          # Relative inhibition to excitation
    K_E = 1000       # Number of excitatory connections
    INH_PROPORTION = 0.25
    V_EXT = 30e-3    # External input potential (V)

    # --- Step 2: Analyze Network State ---
    K_I = K_E * INH_PROPORTION

    # In a balanced network, recurrent excitation and inhibition cancel each other out on average.
    # We check if g * K_I * J is equal to K_E * J.
    # g * K_I = 4.0 * (1000 * 0.25) = 4.0 * 250 = 1000.
    # K_E = 1000.
    # The network is perfectly balanced.
    # Therefore, the mean membrane potential (mu) is determined solely by the external input.
    mu = V_EXT

    # --- Step 3: Formulate Self-Consistent Equations ---
    
    # Integrand for the Siegert formula for firing rate
    def siegert_integrand(x):
        return np.exp(x**2) * (1 + erf(x))

    # Function to calculate firing rate nu for a given sigma
    def calculate_nu(sigma):
        if sigma <= 1e-9: # Avoid division by zero
            # If there's no noise, use the deterministic formula
            if mu <= V_THETA:
                return 0.0
            else:
                isi = TAU_M * np.log((mu - V_RESET) / (mu - V_THETA))
                return 1.0 / (TAU_REF + isi)
        
        # Integration limits
        y_theta = (V_THETA - mu) / sigma
        y_reset = (V_RESET - mu) / sigma

        # Numerical integration
        integral_val, _ = quad(siegert_integrand, y_reset, y_theta)
        
        # Mean inter-spike interval
        isi_mean = TAU_M * np.sqrt(np.pi) * integral_val

        if TAU_REF + isi_mean <= 0:
            # Physically, the rate cannot be infinite. Cap at max possible rate.
            return 1.0 / TAU_REF

        return 1.0 / (TAU_REF + isi_mean)

    # Function to calculate sigma for a given nu
    def calculate_sigma(nu):
        # Variance of the membrane potential due to synaptic bombardment
        sigma_sq = TAU_M * (K_E + G**2 * K_I) * (J**2) * nu
        return np.sqrt(max(0, sigma_sq)) # Ensure non-negative

    # --- Step 4: Solve Numerically ---
    
    # Initial guess for nu using the noise-free case
    nu_current = calculate_nu(0)

    # Iteratively solve for a self-consistent solution
    for i in range(100):
        sigma_current = calculate_sigma(nu_current)
        nu_next_calc = calculate_nu(sigma_current)
        
        # Use a damped update to ensure convergence
        nu_next = 0.9 * nu_current + 0.1 * nu_next_calc

        if abs(nu_next - nu_current) < 1e-4:
            nu_current = nu_next
            break
        nu_current = nu_next
        
    final_nu = nu_current
    final_sigma = calculate_sigma(final_nu)
    
    # --- Step 5: Output the Result ---
    
    print("--- Final Result ---")
    print(f"The converged firing rate is {final_nu:.2f} Hz.")
    print(f"The corresponding standard deviation of the potential is {final_sigma * 1000:.2f} mV.")
    
    # The problem asks to output the numbers in the final equation.
    # The two self-consistent equations are:
    # 1. ν = 1 / (τ_ref + τ * sqrt(π) * ∫[y_reset, y_theta] exp(x²)*(1+erf(x)) dx)
    # 2. σ² = τ * ν * J² * (K_E + g² * K_I)

    print("\n--- Final Equation with Numerical Values ---")

    print("\nEquation 1: Firing Rate (ν)")
    print(f"ν = 1 / (τ_ref + τ * sqrt(π) * Integral)")
    integral_val, _ = quad(siegert_integrand, (V_RESET - mu)/final_sigma, (V_THETA - mu)/final_sigma)
    print(f"{final_nu:.2f} Hz = 1 / ({TAU_REF:.3f} s + {TAU_M:.3f} s * {np.sqrt(np.pi):.3f} * {integral_val:.3f})")
    
    print("\nEquation 2: Potential Fluctuation (σ)")
    print(f"σ² = τ * ν * J² * (K_E + g² * K_I)")
    print(f"({final_sigma:.6f} V)² = {TAU_M:.3f} s * {final_nu:.2f} Hz * ({J:.5f} V)² * ({K_E} + {G**2} * {K_I})")

    final_rate_int = int(round(final_nu))
    print(f"\nThe firing rate of a typical neuron, rounded to the nearest integer, is: {final_rate_int} Hz")
    
    # Returning the final answer in the requested format
    print(f"\n<<<{final_rate_int}>>>")


solve_firing_rate()
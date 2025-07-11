import numpy as np
from scipy.special import erf
from scipy.integrate import quad
from scipy.optimize import fsolve

def solve_and_print_firing_rate():
    """
    Calculates the self-consistent firing rate of a neuron in a random network
    and prints the detailed breakdown of the calculation.
    """
    # Step 1: Define constants from the problem description in base SI units (V, s, Hz).
    tau = 20e-3        # Membrane time constant (s)
    J_E = 0.1e-3       # Excitatory synaptic efficacy (V)
    V_reset = 10e-3    # Reset potential (V)
    V_th = 20e-3       # Threshold potential (V)
    tau_ref = 2e-3     # Refractory period (s)
    g = 4.0            # Relative strength of inhibition to excitation
    K_E = 1000.0       # Number of excitatory connections
    K_I = 0.25 * K_E   # Number of inhibitory connections
    mu_ext = 30e-3     # External input voltage (V)

    # Step 2: Calculate derived parameters.
    J_I = g * J_E      # Inhibitory synaptic efficacy

    # Calculate the mean input potential, mu.
    # The mean recurrent input is proportional to (K_E * J_E - K_I * J_I).
    # (1000 * 0.1e-3) - (250 * 0.4e-3) = 0.1 - 0.1 = 0.
    # The network is perfectly balanced, so the mean is just the external input.
    mu = mu_ext

    # Calculate the variance of the membrane potential, sigma^2.
    # It is linearly dependent on the unknown firing rate nu.
    # sigma^2 = nu * C, where C is the coefficient.
    sigma_squared_coeff = (tau / 2.0) * (K_E * J_E**2 + K_I * J_I**2)

    # Step 3: Define the self-consistent equation to solve for nu.
    # The equation is nu = f(nu), which we rewrite as g(nu) = nu - f(nu) = 0 for the solver.
    def firing_rate_equation(nu):
        # Handle the case of zero firing rate to avoid division by zero.
        if nu <= 1e-9: # Use a small epsilon for stability
            # If rate is zero, noise is zero. Rate is determined by the deterministic formula.
            if mu <= V_th:
                return -0.0 # g(0) = 0 - f(0) = 0 if subthreshold
            T_mean_det = tau_ref + tau * np.log((mu - V_reset) / (mu - V_th))
            nu_det = 1.0 / T_mean_det
            return 0.0 - nu_det

        # Calculate sigma for the current estimate of nu
        sigma = np.sqrt(sigma_squared_coeff * nu)

        # Limits for the Siegert formula integral
        y_th = (V_th - mu) / sigma
        y_reset = (V_reset - mu) / sigma

        # Define the integrand for the formula
        def integrand(u):
            return np.exp(u**2) * (1.0 + erf(u))

        # Numerically compute the integral
        try:
            integral_val, _ = quad(integrand, y_reset, y_th)
        except Exception:
            return np.inf # Return a large number if integration fails

        # Calculate the new firing rate based on the integral
        T_mean = tau_ref + tau * np.sqrt(np.pi) * integral_val
        if T_mean <= 0: return np.inf
        
        nu_new = 1.0 / T_mean

        # Return the difference for the root-finder
        return nu - nu_new

    # Step 4: Solve the equation numerically.
    # We use the deterministic rate as a good initial guess.
    T_det = tau_ref + tau * np.log((mu - V_reset) / (mu - V_th))
    nu_initial_guess = 1.0 / T_det
    
    # Use fsolve to find the root
    nu_final = fsolve(firing_rate_equation, nu_initial_guess)[0]

    # Step 5: Print the detailed results as requested.
    # Recalculate components with the final nu for printing.
    sigma_final = np.sqrt(sigma_squared_coeff * nu_final)
    y_th_final = (V_th - mu) / sigma_final
    y_reset_final = (V_reset - mu) / sigma_final
    integral_val_final, _ = quad(lambda u: np.exp(u**2) * (1.0 + erf(u)), y_reset_final, y_th_final)
    T_mean_final = tau_ref + tau * np.sqrt(np.pi) * integral_val_final

    print("The firing rate ν is found by solving the self-consistent equation:")
    print("ν = 1 / (τ_ref + τ * sqrt(π) * ∫[y_reset, y_th] exp(u²) * (1 + erf(u)) du)\n")
    print(f"The converged firing rate is ν = {nu_final:.2f} Hz.\n")
    print("For this rate, the parameters in the final equation are:")
    print("-" * 50)
    print(f"Membrane Time Constant (τ):     {tau*1000:.0f} ms")
    print(f"Refractory Period (τ_ref):        {tau_ref*1000:.0f} ms")
    print(f"Voltage Threshold (V_th):       {V_th*1000:.0f} mV")
    print(f"Voltage Reset (V_reset):        {V_reset*1000:.0f} mV")
    print(f"Mean Input (μ):                 {mu*1000:.0f} mV")
    print(f"Fluctuation Strength (σ):       {sigma_final*1000:.3f} mV (from ν = {nu_final:.2f} Hz)")
    print(f"Lower Integration Limit (y_reset):{y_reset_final:.4f}")
    print(f"Upper Integration Limit (y_th):   {y_th_final:.4f}")
    print(f"Value of Integral:                {integral_val_final:.4f}")
    print(f"Mean Inter-Spike Interval (T):  {T_mean_final*1000:.2f} ms")
    print("-" * 50)
    print(f"Final check: 1 / T = 1 / {T_mean_final:.4f} s = {1/T_mean_final:.2f} Hz, which matches the converged ν.")
    
    # The problem asks for the answer as an integer.
    final_answer = int(round(nu_final))
    print(f"\nRounding the result gives a final firing rate of {final_answer} Hz.")
    print(f"\n<<<{final_answer}>>>")

solve_and_print_firing_rate()
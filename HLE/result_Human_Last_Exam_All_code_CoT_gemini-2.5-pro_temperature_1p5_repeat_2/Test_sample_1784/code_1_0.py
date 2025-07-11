import numpy as np
from scipy.integrate import quad
from scipy.special import erfc
from scipy.optimize import fsolve

def solve_firing_rate():
    """
    Calculates the firing rate of a neuron in a balanced network.
    """
    # --- Step 1: Define Parameters in SI units ---
    tau_m = 0.020   # Membrane time constant (20 ms)
    V_th = 0.020    # Voltage threshold (20 mV)
    V_reset = 0.010 # Voltage reset (10 mV)
    tau_ref = 0.002 # Refractory period (2 ms)
    J_E = 0.0001    # Excitatory synaptic efficacy (0.1 mV)
    g = 4.0         # Relative inhibition to excitation
    J_I = g * J_E   # Inhibitory synaptic efficacy (0.4 mV)
    K_E = 1000      # Number of excitatory connections
    K_I = 250       # Number of inhibitory connections (0.25 * K_E)
    mu_ext = 0.030  # External input drive (30 mV)

    # --- Step 2: Calculate Mean and Variance of Input ---
    # The network is in the balanced state because K_E*J_E (1000*0.1mV=100mV)
    # equals K_I*J_I (250*0.4mV=100mV).
    # Therefore, the mean recurrent input is zero.
    mu = mu_ext

    # The variance of the membrane potential is sigma_V^2 = C * nu, where nu is the firing rate.
    # C = (tau_m / 2.0) * (K_E * J_E^2 + K_I * J_I^2)
    C = (tau_m / 2.0) * (K_E * J_E**2 + K_I * J_I**2)

    # --- Step 3: Formulate the Self-Consistent Equation ---
    # The integrand for the mean first passage time calculation
    def integrand(u):
        return np.exp(u**2) * erfc(-u)

    # This function represents the self-consistent equation: f(nu) - nu = 0
    # where f(nu) is the calculated output rate for a given network rate nu.
    def rate_equation(nu):
        # Firing rate must be positive
        if nu <= 0:
            return 1.0 # Return a non-zero value to guide the solver

        # Calculate sigma_V from the input rate nu
        sigma_V_sq = C * nu
        if sigma_V_sq <= 0:
             return 1.0
        sigma_V = np.sqrt(sigma_V_sq)

        # Calculate integral limits y_reset and y_th
        y_th = (V_th - mu) / sigma_V
        y_reset = (V_reset - mu) / sigma_V

        # Calculate the mean first passage time T
        # Since mu > V_th, the neuron is in the supra-threshold regime.
        integral_val, integral_err = quad(integrand, y_reset, y_th)
        T = tau_m * np.sqrt(np.pi) * integral_val
        
        # T should be positive. If not, it's a numerical artifact.
        if T <= 0:
            # Fallback to deterministic formula for robustness
            if mu > V_th:
                T = tau_m * np.log((mu - V_reset) / (mu - V_th))
            else:
                 # Should not happen in this problem
                return 1.0

        # Calculate the resulting output rate
        nu_out = 1.0 / (tau_ref + T)

        return nu_out - nu

    # --- Step 4: Solve the Equation ---
    # Use the deterministic rate (no noise) as a good initial guess
    if mu > V_th:
        T_det = tau_m * np.log((mu - V_reset) / (mu - V_th))
        nu_guess = 1.0 / (tau_ref + T_det)
    else:
        nu_guess = 10.0 # A small positive guess if subthreshold

    # Use a numerical solver to find the root of the rate_equation
    firing_rate_solution = fsolve(rate_equation, x0=nu_guess)
    final_rate_hz = firing_rate_solution[0]

    # --- Step 5: Present the Final Answer ---
    print(f"The self-consistent firing rate is {final_rate_hz:.2f} Hz.")
    final_rate_int = int(round(final_rate_hz))
    print(f"The firing rate of a typical neuron is {final_rate_int} Hz.")
    print("\n--- Final Equation Details ---")
    print("The firing rate ν is the solution to the self-consistent equation: ν = 1 / (τ_ref + T(ν))")
    
    # Recalculate values at the solution point for printing
    sigma_V_final = np.sqrt(C * final_rate_hz)
    y_th_final = (V_th - mu) / sigma_V_final
    y_reset_final = (V_reset - mu) / sigma_V_final
    integral_final, _ = quad(integrand, y_reset_final, y_th_final)
    T_final = tau_m * np.sqrt(np.pi) * integral_final
    
    print(f"Solved for ν = {final_rate_hz:.2f} Hz, where the following values hold:")
    print(f"  τ_ref (refractory period) = {tau_ref * 1000:.1f} ms")
    print(f"  T(ν) (mean time-to-threshold) = {T_final * 1000:.2f} ms")
    print("T(ν) is calculated using T = τ_m * sqrt(π) * ∫[y_reset, y_th] g(u)du, with:")
    print(f"  τ_m (membrane time constant) = {tau_m * 1000:.1f} ms")
    print(f"  μ (mean input) = {mu * 1000:.1f} mV")
    print(f"  σ_V (input fluctuation std dev) = {sigma_V_final * 1000:.2f} mV")
    print(f"  V_th (threshold) = {V_th * 1000:.1f} mV")
    print(f"  V_reset (reset) = {V_reset * 1000:.1f} mV")
    print(f"  y_reset ((V_reset - μ)/σ_V) = {y_reset_final:.2f}")
    print(f"  y_th ((V_th - μ)/σ_V) = {y_th_final:.2f}")
    print(f"Check: 1 / ({tau_ref:.3f} s + {T_final:.4f} s) = {1/(tau_ref + T_final):.2f} Hz, which matches the solution.")

if __name__ == '__main__':
    solve_firing_rate()
<<<62>>>
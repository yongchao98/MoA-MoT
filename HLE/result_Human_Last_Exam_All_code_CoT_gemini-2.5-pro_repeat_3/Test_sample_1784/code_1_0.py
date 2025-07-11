import numpy as np
from scipy import integrate, special, optimize

def solve_firing_rate():
    """
    Calculates the firing rate of a neuron in a balanced network.
    """
    # 1. Define network parameters in SI units
    tau_m = 20e-3      # Membrane time constant (s)
    J = 0.1e-3         # Excitatory synaptic efficacy (V)
    V_reset = 10e-3    # Reset voltage (V)
    V_th = 20e-3       # Voltage threshold (V)
    tau_ref = 2e-3     # Refractory period (s)
    g = 4.0            # Relative strength of inhibition to excitation
    K_E = 1000         # Number of excitatory connections
    prop_I = 0.25      # Proportion of inhibitory connections relative to K_E
    V_ext = 30e-3      # External input potential (V)

    # Derived parameters
    K_I = prop_I * K_E
    J_E = J
    J_I = g * J_E

    # 2. Calculate mean (mu) and variance factor for the membrane potential
    # In a balanced network, the mean recurrent input is zero.
    # mu = V_ext + tau_m * nu * (K_E * J_E - K_I * J_I)
    # The term (K_E * J_E - K_I * J_I) = (1000 * 0.1e-3 - 250 * 0.4e-3) = 0.1 - 0.1 = 0.
    mu = V_ext

    # The variance of the potential is sigma^2 = factor * nu
    # factor = (tau_m / 2) * (K_E * J_E^2 + K_I * J_I^2)
    variance_factor = (tau_m / 2) * (K_E * J_E**2 + K_I * J_I**2)

    # 3. Define the self-consistency equation: F(nu) = nu - calculated_nu(nu) = 0
    def self_consistency_equation(nu):
        """
        Represents the equation nu = f(nu), which we need to solve.
        Returns the difference nu - f(nu).
        """
        # nu is the firing rate in Hz. It must be positive.
        if nu <= 0:
            return -1

        # Calculate sigma for the current guess of nu
        sigma = np.sqrt(variance_factor * nu)

        # The firing rate formula is unstable for zero noise, use the deterministic limit.
        if sigma < 1e-7:
            if mu > V_th:
                # Deterministic firing rate for supra-threshold input
                T_charge = tau_m * np.log((mu - V_reset) / (mu - V_th))
                calculated_nu = 1.0 / (tau_ref + T_charge)
                return nu - calculated_nu
            else:
                return nu # Firing rate is zero for sub-threshold input without noise

        # Define the integral from the theoretical firing rate formula
        y_th = (V_th - mu) / sigma
        y_reset = (V_reset - mu) / sigma
        
        integrand = lambda y: np.exp(y**2) * (1 + special.erf(y))
        
        try:
            integral_val, _ = integrate.quad(integrand, y_reset, y_th)
        except Exception:
            # Return a large number if integration fails to guide the solver
            return np.inf

        # Calculate the mean first passage time from V_reset to V_th
        T_passage = tau_m * np.sqrt(np.pi) * integral_val

        # If T_passage is negative (can happen for large mu), the neuron fires instantly.
        if (tau_ref + T_passage) <= 0:
            return np.inf # Effectively, rate is infinite, should not happen here.

        calculated_nu = 1.0 / (tau_ref + T_passage)
        
        return nu - calculated_nu

    # 4. Solve the equation numerically
    # Provide an initial guess using the deterministic formula (ignoring noise)
    T_charge_det = tau_m * np.log((mu - V_reset) / (mu - V_th))
    nu_guess = 1.0 / (tau_ref + T_charge_det)

    # Use a robust root-finding algorithm
    try:
        solution = optimize.fsolve(self_consistency_equation, x0=nu_guess)
        final_rate_hz = solution[0]
    except Exception as e:
        print(f"Solver failed: {e}. The result might not be accurate.")
        final_rate_hz = 0

    # 5. Output the results
    # Recalculate final values for printing
    sigma_final = np.sqrt(variance_factor * final_rate_hz)
    y_th_final = (V_th - mu) / sigma_final
    y_reset_final = (V_reset - mu) / sigma_final
    integrand_final = lambda y: np.exp(y**2) * (1 + special.erf(y))
    integral_val_final, _ = integrate.quad(integrand_final, y_reset_final, y_th_final)
    T_passage_final = tau_m * np.sqrt(np.pi) * integral_val_final

    print("The firing rate is determined by the equation: Rate = 1 / (Refractory_Period + Passage_Time)")
    print("Using the solved rate, the values are:")
    print(f"Refractory_Period = {tau_ref:.3f} s")
    print(f"Passage_Time = {T_passage_final:.3f} s")
    print("\nFinal equation:")
    print(f"{int(round(final_rate_hz))} = 1 / ({tau_ref:.3f} + {T_passage_final:.3f})")
    print(f"\nThe firing rate of a typical neuron is {int(round(final_rate_hz))} Hz.")

    return int(round(final_rate_hz))

# Run the solver and store the final answer
final_answer = solve_firing_rate()
print(f'<<<{final_answer}>>>')
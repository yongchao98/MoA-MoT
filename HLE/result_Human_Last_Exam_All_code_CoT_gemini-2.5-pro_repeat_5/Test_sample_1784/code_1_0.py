import numpy as np
from scipy.special import erf
from scipy.integrate import quad
from scipy.optimize import root_scalar

def solve_firing_rate():
    """
    Calculates the firing rate of a neuron in a balanced random network.
    """
    # Step 1: Define network and neuron parameters in SI units
    tau_m = 20e-3      # Membrane time constant (s)
    J_efficacy = 0.1e-3  # Synaptic efficacy (V)
    V_reset = 10e-3    # Reset voltage (V)
    V_th = 20e-3       # Threshold voltage (V)
    tau_ref = 2e-3     # Refractory period (s)
    g = 4.0            # Relative strength of inhibition to excitation
    K_E = 1000         # Number of excitatory connections
    K_I = 250          # Number of inhibitory connections (0.25 * K_E)
    mu_ext = 30e-3     # External input voltage (V)

    # Step 2: Define derived parameters
    J_E = J_efficacy
    J_I = -g * J_E

    # Step 3: Define mean and standard deviation of the membrane potential.
    # The mean input from the recurrent network is tau_m * r * (K_E * J_E + K_I * J_I)
    # K_E * J_E = 1000 * 0.1e-3 = 0.1 V
    # K_I * J_I = 250 * (-4 * 0.1e-3) = -0.1 V
    # The sum is 0, so the network is in a balanced state.
    # The mean potential is constant and equals the external input.
    mu = mu_ext

    # The variance of the potential depends on the firing rate 'r'.
    # sigma^2 = tau_m * r * (K_E * J_E^2 + K_I * J_I^2)
    sigma_sq_factor = tau_m * (K_E * J_E**2 + K_I * J_I**2)

    def get_sigma(r):
        """Calculates sigma from the firing rate r."""
        # Add a small epsilon to r to prevent sqrt(0) if r is exactly 0
        return np.sqrt(sigma_sq_factor * (r + 1e-12))

    # The integrand in the Siegert formula for the mean first passage time
    def integrand(x):
        return np.exp(x**2) * (1 + erf(x))

    # Step 4: Define the self-consistent equation G(r) = r - f(r) = 0
    def firing_rate_equation(r):
        """
        Represents the self-consistent equation for the firing rate r.
        We want to find the root of this function, i.e., where it equals zero.
        """
        if r < 0:
            return 1.0 # Rates must be non-negative

        sigma = get_sigma(r)
        
        # Integration bounds for the Siegert formula
        u_reset = (V_reset - mu) / sigma
        u_thresh = (V_th - mu) / sigma

        # The Siegert formula can be numerically unstable for large negative inputs.
        # For large negative u, T_passage approximates the noiseless case.
        if u_thresh < -35: # exp(35^2) overflows standard floats
            T_passage = tau_m * np.log((V_reset - mu) / (V_th - mu))
        else:
            integral_val, _ = quad(integrand, u_reset, u_thresh)
            T_passage = tau_m * np.sqrt(np.pi) * integral_val
        
        # Calculate the new firing rate f(r) based on the total time to fire
        total_time = tau_ref + T_passage
        
        # Avoid division by zero if total_time is vanishingly small
        if total_time <= 1e-9:
            r_new = 1.0 / tau_ref
        else:
            r_new = 1.0 / total_time
        
        return r - r_new

    # Step 5: Solve the equation G(r) = 0 numerically
    # The firing rate is bounded by 0 and 1/tau_ref.
    # We use a robust root-finding method within this bracket.
    try:
        sol = root_scalar(firing_rate_equation, bracket=[0.1, 1/tau_ref - 1], method='brentq')
        final_rate = sol.root
    except ValueError:
        print("Error: Could not find a solution for the firing rate in the given bracket.")
        return

    # Step 6: Output the final result and the final equation numbers
    print(f"The calculated firing rate of a typical neuron is {int(round(final_rate))} Hz.")
    print("\nThis result is the solution to the self-consistent equation r = 1 / (τ_ref + T_passage(r)).")
    print("Using the final firing rate to calculate the terms in the equation gives:")
    
    # Recalculate final values for printing
    sigma_final = get_sigma(final_rate)
    u_reset_final = (V_reset - mu) / sigma_final
    u_thresh_final = (V_th - mu) / sigma_final
    integral_val_final, _ = quad(integrand, u_reset_final, u_thresh_final)
    T_passage_final = tau_m * np.sqrt(np.pi) * integral_val_final

    print(f"\nFinal Firing Rate (r): {final_rate:.2f} Hz")
    print(f"Membrane Time Constant (τ_m): {tau_m} s")
    print(f"Refractory Period (τ_ref): {tau_ref} s")
    print(f"Mean Input Voltage (μ): {mu} V")
    print(f"Std. Dev. of Input Voltage (σ): {sigma_final:.6f} V")
    print(f"Mean First Passage Time (T_passage): {T_passage_final:.6f} s")
    
    print("\nFinal Equation:")
    print(f"{final_rate:.2f} = 1 / ({tau_ref} + {T_passage_final:.6f})")
    
    # Store final integer answer for automated checking
    global final_answer
    final_answer = int(round(final_rate))


if __name__ == '__main__':
    final_answer = None
    solve_firing_rate()
    # The final answer format is specified by the user.
    # print(f"<<<{final_answer}>>>")
import numpy as np
from scipy.integrate import quad
from scipy.special import erf

def solve_firing_rate():
    """
    This function calculates the firing rate of a neuron in a balanced recurrent network.
    """
    # 1. Define network and neuron parameters in SI units.
    tau_m = 20e-3      # Membrane time constant (s)
    J_E = 0.1e-3       # Excitatory synaptic efficacy (V)
    V_reset = 10e-3    # Reset potential (V)
    V_th = 20e-3       # Threshold potential (V)
    tau_ref = 2e-3     # Refractory period (s)
    g = 4.0            # Relative strength of inhibition
    K_E = 1000         # Number of excitatory connections
    K_I = 250          # Number of inhibitory connections
    V_ext = 30e-3      # External input voltage (V)

    J_I = g * J_E      # Inhibitory synaptic efficacy (V)

    # 2. Calculate mean and standard deviation of the membrane potential.
    # The network is in a balanced state because the mean recurrent input is zero:
    # K_E * J_E = 1000 * 0.1e-3 = 0.1 V
    # K_I * J_I = 250 * (4 * 0.1e-3) = 0.1 V
    # Thus, the mean potential is determined only by the external input.
    mu_V = V_ext

    # The variance of the potential fluctuations depends on the firing rate 'nu'.
    # σ_V^2 = ν * [τ_m * (K_E*J_E^2 + K_I*J_I^2)]
    var_coeff = tau_m * (K_E * J_E**2 + K_I * J_I**2)

    # 3. & 4. Solve the self-consistent equation for the firing rate 'nu'.
    # This function defines the integrand in Siegert's formula for the mean first passage time.
    def siegert_integrand(x):
        return np.exp(x**2) * (1 + erf(x))

    # This function computes the firing rate 'nu_out' for a given input rate 'nu_in'.
    # We will iterate this function to find the fixed point where nu_out = nu_in.
    def calculate_nu_from_input(nu_in):
        # Handle the case of zero noise (or zero input rate).
        if nu_in <= 0:
            if mu_V > V_th:
                # If mean drive is above threshold, firing is deterministic.
                t_isi_noiseless = tau_m * np.log((mu_V - V_reset) / (mu_V - V_th))
                return 1.0 / (tau_ref + t_isi_noiseless)
            else:
                return 0.0

        sigma_V = np.sqrt(var_coeff * nu_in)
        if sigma_V < 1e-9: # Effectively zero noise
             return calculate_nu_from_input(0)

        # Normalized integration limits for the Siegert formula
        y_th = (V_th - mu_V) / sigma_V
        y_reset = (V_reset - mu_V) / sigma_V

        # Numerically compute the integral
        integral_val, _ = quad(siegert_integrand, y_reset, y_th)

        # Calculate the effective mean inter-spike interval (ISI)
        t_isi_eff = tau_m * np.sqrt(np.pi) * integral_val

        # The total period is ISI + refractory period. Rate is its inverse.
        return 1.0 / (tau_ref + t_isi_eff)

    # Use the noiseless rate as an initial guess for the iteration.
    nu = calculate_nu_from_input(0)

    # Iteratively solve for nu until convergence.
    for _ in range(100):
        nu_new = calculate_nu_from_input(nu)
        if np.abs(nu_new - nu) < 1e-4:
            nu = nu_new
            break
        nu = nu_new

    # 5. Present the final result.
    final_rate = nu
    
    # Recalculate the final ISI with the converged rate for printing.
    final_t_isi = (1.0 / final_rate) - tau_ref
    
    print(f"The firing rate of a typical neuron can be calculated as the inverse of the total spike period.")
    print(f"Total Period = Refractory Period + Mean Inter-Spike Interval (ISI)")
    print(f"\nFinal Equation:")
    print(f"{int(round(final_rate))} Hz = 1 / ({tau_ref:.3f} s + {final_t_isi:.4f} s)")

solve_firing_rate()
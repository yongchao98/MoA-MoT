import numpy as np
from scipy.special import erf
from scipy.integrate import quad
from scipy.optimize import fsolve

def solve_firing_rate():
    """
    Calculates the firing rate of a neuron in a randomly connected network.
    """
    # --- 1. Define Network and Neuron Parameters ---
    # Using SI units (seconds, Volts) for internal calculations
    tau_m = 20e-3      # Membrane time constant (s)
    J = 0.1e-3         # Synaptic efficacy (V)
    V_reset = 10e-3    # Voltage reset (V)
    V_th = 20e-3       # Voltage threshold (V)
    tau_ref = 2e-3     # Refractory period (s)
    g = 4.0            # Relative strength of inhibition to excitation
    K_E = 1000         # Number of excitatory connections
    K_I = 250          # Number of inhibitory connections (0.25 * K_E)
    mu_ext = 30e-3     # External input mean potential (V)

    # Synaptic weights
    J_E = J
    J_I = -g * J_E

    # --- 2. Calculate Mean and Variance of Membrane Potential ---
    # The mean recurrent input is K_E*J_E + K_I*J_I = 1000*0.1e-3 + 250*(-0.4e-3) = 0.1 - 0.1 = 0.
    # The network is balanced, so the mean potential is driven only by external input.
    mu_V = mu_ext

    # The variance of the potential (sigma_V^2) is proportional to the firing rate nu.
    # sigma_V^2 = C_sigma_sq * nu
    # The factor of 0.5 comes from the integral of the synaptic filter correlation function.
    C_sigma_sq = (tau_m / 2) * (K_E * J_E**2 + K_I * J_I**2)

    # --- 3. Define the Self-Consistent Equation for Firing Rate ---
    # We need to solve g(nu) = f(nu) - nu = 0, where f(nu) is the theoretical rate.
    def firing_rate_equation(nu):
        # Avoid division by zero or math domain errors if nu is non-positive
        if nu <= 0:
            # If there's no firing, there's no variance (noise-free case)
            # In the suprathreshold regime (mu_V > V_th), the rate is non-zero.
            if mu_V > V_th:
                T_isi = tau_m * np.log((mu_V - V_reset) / (mu_V - V_th))
                return 1.0 / (tau_ref + T_isi)
            else:
                return 0.0

        # Calculate variance for the given firing rate nu
        sigma_V_sq = C_sigma_sq * nu
        sigma_V = np.sqrt(sigma_V_sq)

        # The firing rate is given by the Siegert formula for the mean first passage time.
        # nu = 1 / (tau_ref + T_isi)
        # T_isi = tau_m * sqrt(pi) * integral
        
        # Define the integrand for the T_isi calculation
        # This is the core of the Siegert formula
        integrand = lambda x: np.exp(x**2) * (1 + erf(x))

        # Define the integration limits
        y_th = (V_th - mu_V) / sigma_V
        y_reset = (V_reset - mu_V) / sigma_V

        # Numerically compute the integral
        integral_val, _ = quad(integrand, y_reset, y_th)

        # Calculate the mean inter-spike interval (T_isi)
        T_isi = tau_m * np.sqrt(np.pi) * integral_val
        
        # Return the calculated firing rate
        if (tau_ref + T_isi) <= 0:
            return np.inf # Should not happen with these parameters
        return 1.0 / (tau_ref + T_isi)

    # --- 4. Solve the Equation Numerically ---
    # Define the function to find the root of: f(nu) - nu = 0
    equation_to_solve = lambda nu: firing_rate_equation(nu) - nu

    # Initial guess: use the noise-free rate
    initial_guess = 1.0 / (tau_ref + tau_m * np.log((mu_V - V_reset) / (mu_V - V_th)))
    
    # Solve for the equilibrium firing rate
    nu_final = fsolve(equation_to_solve, x0=initial_guess)[0]

    # --- 5. Print the Results ---
    print(f"The equilibrium firing rate of a typical neuron is {int(round(nu_final))} Hz.\n")
    print("This result is derived from solving the self-consistent equation for the network.")
    print("Below are the final values that satisfy the equation:\n")

    # Recalculate final values for printing
    sigma_V_final = np.sqrt(C_sigma_sq * nu_final)
    y_th_final = (V_th - mu_V) / sigma_V_final
    y_reset_final = (V_reset - mu_V) / sigma_V_final
    integrand = lambda x: np.exp(x**2) * (1 + erf(x))
    integral_final, _ = quad(integrand, y_reset_final, y_th_final)
    T_isi_final = tau_m * np.sqrt(np.pi) * integral_final
    
    print("Final Equation: nu = 1 / (tau_ref + T_isi)")
    print("-------------------------------------------------")
    print(f"Converged Firing Rate (nu): {nu_final:.2f} Hz")
    print(f"Refractory Period (tau_ref): {tau_ref * 1000:.2f} ms")
    print(f"Mean Inter-Spike Interval (T_isi): {T_isi_final * 1000:.2f} ms")
    print(f"Total Period (tau_ref + T_isi): {(tau_ref + T_isi_final) * 1000:.2f} ms")
    print(f"Calculated Rate (1 / Total Period): {1/(tau_ref + T_isi_final):.2f} Hz\n")

    print("T_isi Calculation: T_isi = tau_m * sqrt(pi) * integral")
    print("-------------------------------------------------")
    print(f"Membrane Time Constant (tau_m): {tau_m * 1000:.2f} ms")
    print(f"Integral Value: {integral_final:.4f}\n")

    print("Integral Bounds: y = (V - mu_V) / sigma_V")
    print("-------------------------------------------------")
    print(f"Mean Potential (mu_V): {mu_V * 1000:.2f} mV")
    print(f"Std Dev of Potential (sigma_V): {sigma_V_final * 1000:.2f} mV")
    print(f"Threshold Voltage (V_th): {V_th * 1000:.2f} mV")
    print(f"Reset Voltage (V_reset): {V_reset * 1000:.2f} mV")
    print(f"Upper Bound (y_th): {y_th_final:.2f}")
    print(f"Lower Bound (y_reset): {y_reset_final:.2f}")
    
    return int(round(nu_final))

if __name__ == '__main__':
    final_rate = solve_firing_rate()
    # The final answer is submitted in the special format below
    # print(f"\n<<< {final_rate} >>>")
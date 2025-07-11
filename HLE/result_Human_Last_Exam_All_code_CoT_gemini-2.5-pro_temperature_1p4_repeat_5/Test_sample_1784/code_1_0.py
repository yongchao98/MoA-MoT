import math

def solve_firing_rate():
    """
    Calculates the firing rate of a neuron in a balanced network using a self-consistent mean-field theory.
    """
    # 1. Define parameters in standard units (Volts, Seconds)
    tau_m = 0.020       # Membrane time constant (s)
    J = 0.0001          # Synaptic efficacy (V)
    V_reset = 0.010     # Reset potential (V)
    V_th = 0.020        # Threshold potential (V)
    tau_ref = 0.002     # Refractory period (s)
    g = 4.0             # Relative strength of inhibition
    K_E = 1000          # Number of excitatory connections
    
    # Inhibitory connections are 0.25 of excitatory connections
    K_I = 0.25 * K_E    
    
    # External input sets the mean potential in a balanced network
    mu_v = 0.030        # Mean membrane potential (V)

    # 2. Check for the balanced condition
    # This should be zero, simplifying the mean potential calculation
    balance_check = K_E - g * K_I

    # This term is used in the variance calculation
    variance_coeff = K_E + g**2 * K_I

    # Helper function to calculate the standard deviation of Vm from the firing rate nu
    def get_sigma_v(nu):
        # sigma_V^2 = (tau_m/2) * J^2 * nu * (K_E + g^2 * K_I)
        if nu <= 0:
            return 1e-9 # Return a very small number to avoid errors, cannot be 0
        sigma_v_sq = (tau_m / 2.0) * (J**2) * nu * variance_coeff
        return math.sqrt(sigma_v_sq)

    # Helper function for the integrand in the Siegert formula
    def integrand(x):
        return math.exp(x**2) * (1 + math.erf(x))

    # Helper function for numerical integration (Trapezoidal Rule)
    def integrate_siegert(y_lo, y_hi, n_steps=1000):
        if y_lo >= y_hi:
            return 0.0
        h = (y_hi - y_lo) / n_steps
        integral = 0.5 * (integrand(y_lo) + integrand(y_hi))
        for i in range(1, n_steps):
            integral += integrand(y_lo + i * h)
        return integral * h

    # 3. Solve for nu using fixed-point iteration
    # Initial guess for nu (e.g., deterministic rate when mean is supra-threshold)
    try:
        t_isi_det = tau_m * math.log((mu_v - V_reset) / (mu_v - V_th))
        nu = 1.0 / (tau_ref + t_isi_det)
    except ValueError:
        nu = 50.0 # Fallback guess

    # Iterate to find the self-consistent solution
    for _ in range(100):
        sigma_v = get_sigma_v(nu)

        # Lower and upper bounds for the integral
        y_hi = (V_th - mu_v) / sigma_v
        y_lo = (V_reset - mu_v) / sigma_v

        integral_val = integrate_siegert(y_lo, y_hi)
        
        # This is the mean time to fire from reset, based on the Siegert solution
        tau_s = tau_m * math.sqrt(math.pi) * integral_val

        if (tau_ref + tau_s) <= 1e-9: # Avoid division by zero
             # This state indicates an extremely high firing rate
             nu_new = 1000.0
        else:
             nu_new = 1.0 / (tau_ref + tau_s)

        # Check for convergence
        if abs(nu_new - nu) < 1e-4:
            nu = nu_new
            break
        nu = nu_new
    
    # 4. Output the final calculation and result
    final_rate = nu
    final_rate_int = int(round(final_rate))
    
    sigma_v_final = get_sigma_v(final_rate)
    y_hi_final = (V_th - mu_v) / sigma_v_final
    y_lo_final = (V_reset - mu_v) / sigma_v_final
    integral_val_final = integrate_siegert(y_lo_final, y_hi_final)
    tau_s_final = tau_m * math.sqrt(math.pi) * integral_val_final

    print("Final firing rate equation: ν = 1 / (τ_ref + τ_s)")
    print("Where τ_s = τ_m * √π * ∫[y_lo, y_hi] exp(x²)(1+erf(x)) dx\n")
    print(f"Converged Firing Rate (ν): {final_rate:.2f} Hz")
    print("---------------------------------")
    print("Final Calculation Values:")
    print(f"  τ_ref = {tau_ref:.3f} s")
    print(f"  τ_s = {tau_s_final:.5f} s")
    print(f"    τ_m = {tau_m:.3f} s")
    print(f"    Integral Value = {integral_val_final:.4f}")
    print(f"    Integration bounds: y_lo={y_lo_final:.2f}, y_hi={y_hi_final:.2f}")
    print("---------------------------------")
    print(f"{final_rate:.2f} = 1 / ({tau_ref:.3f} + {tau_s_final:.5f})")

    print(f"\nThe final firing rate is approximately {final_rate_int} Hz.")
    print(f"<<<{final_rate_int}>>>")

solve_firing_rate()
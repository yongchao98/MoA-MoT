import numpy as np
from scipy.integrate import quad
from scipy.special import erf
from scipy.optimize import fsolve

# --- Define neuron and network parameters ---
tau_m_ms = 20.0  # Membrane time constant (ms)
V_th_mv = 20.0   # Voltage threshold (mV)
V_reset_mv = 10.0 # Voltage reset (mV)
tau_ref_ms = 2.0  # Refractory period (ms)
g = 4.0           # Relative inhibition to excitation
K_E = 1000        # Number of excitatory connections
J_E_mv = 0.1      # Excitatory synaptic efficacy (mV)
mu_ext_mv = 30.0  # External input drive (mV)
K_I_ratio = 0.25  # Proportion of inhibitory to excitatory connections

# Convert to SI units for calculations
tau_m_s = tau_m_ms / 1000.0
tau_ref_s = tau_ref_ms / 1000.0

# --- Calculate µ and the factor for σ ---
# In this balanced network, the mean recurrent input is zero.
K_I = K_E * K_I_ratio
J_I_mv = g * J_E_mv
mean_recurrent_factor = K_E * J_E_mv - K_I * J_I_mv

# Total mean potential µ is constant
mu_V_mv = mu_ext_mv

# Calculate the factor for the variance of the potential
# sigma_V^2 = (tau_m * (K_E*J_E^2 + K_I*J_I^2)) * r = sigma_factor_sq * r
# Units: s * (mV^2) * Hz -> mV^2
sigma_factor_sq = tau_m_s * (K_E * J_E_mv**2 + K_I * J_I_mv**2)

# --- Define the self-consistency equation ---
# The equation to solve is f(r) = r - theoretical_rate(r) = 0
def firing_rate_equation(r):
    """
    Defines the equation f(r) = r - F(µ, σ(r)).
    The root of this equation is the self-consistent firing rate.
    'r' is assumed to be a single-element array as required by fsolve.
    """
    rate_hz = r[0]

    # Avoid math errors for non-physical rates
    if rate_hz <= 0:
        return 1.0

    # Calculate standard deviation σ from the firing rate r
    sigma_V_sq = sigma_factor_sq * rate_hz
    sigma_V_mv = np.sqrt(sigma_V_sq)
    
    # Calculate integration bounds for the firing rate formula
    y_th = (V_th_mv - mu_V_mv) / sigma_V_mv
    y_reset = (V_reset_mv - mu_V_mv) / sigma_V_mv

    # Define the integrand for the Siegert formula
    integrand = lambda y: np.exp(y**2) * (1 + erf(y))
    
    # Perform the numerical integration
    try:
        integral_val, _ = quad(integrand, y_reset, y_th)
    except Exception as e:
        print(f"Integration failed: {e}")
        return np.inf

    # Calculate the theoretical mean time to first passage
    t_fp_s = tau_m_s * np.sqrt(np.pi) * integral_val
    
    # Calculate theoretical firing rate in Hz
    # Add a small epsilon to the denominator to prevent division by zero
    # This might happen if t_fp becomes -tau_ref, which is unlikely here.
    denominator = tau_ref_s + t_fp_s
    if denominator <= 0:
        return np.inf
        
    theoretical_rate_hz = 1.0 / denominator
    
    return rate_hz - theoretical_rate_hz

# --- Solve the equation ---
# Initial guess for the firing rate (e.g., from a noiseless model)
T_isi_noiseless = tau_m_s * np.log((mu_V_mv - V_reset_mv) / (mu_V_mv - V_th_mv))
initial_guess = [1.0 / (tau_ref_s + T_isi_noiseless)]

# Use fsolve to find the root of the equation
solution_rate, info, ier, mesg = fsolve(firing_rate_equation, initial_guess, full_output=True)

if ier == 1:
    final_rate_hz = solution_rate[0]
    final_rate_int = int(round(final_rate_hz))

    print("--- Model Parameters and Derived Values ---")
    print(f"Mean potential µ = {mu_V_mv:.1f} mV")
    print(f"Variance of potential σ² = {sigma_factor_sq:.2f} * r (mV²)")
    print("")
    print("--- Self-Consistency Solution ---")
    print(f"Solved the equation: r = 1 / (τ_ref + τ_m * √π * ∫[exp(y²) * (1 + erf(y))] dy)")
    print("The solved firing rate and the numbers used in the final equation are:")
    print(f"  Firing Rate (r): {final_rate_hz:.2f} Hz, rounded to {final_rate_int} Hz")
    
    # Recalculate values at the solution for printing
    sigma_final = np.sqrt(sigma_factor_sq * final_rate_hz)
    y_th_final = (V_th_mv - mu_V_mv) / sigma_final
    y_reset_final = (V_reset_mv - mu_V_mv) / sigma_final
    
    print("\n--- Final Equation Parameters ---")
    print(f"  Voltage Threshold (V_th): {V_th_mv:.1f} mV")
    print(f"  Reset Voltage (V_reset): {V_reset_mv:.1f} mV")
    print(f"  Mean Potential (µ): {mu_V_mv:.1f} mV")
    print(f"  Std. Dev. of Potential (σ): {sigma_final:.2f} mV")
    print(f"  Membrane Time Constant (τ_m): {tau_m_s:.3f} s")
    print(f"  Refractory Period (τ_ref): {tau_ref_s:.3f} s")
    print(f"  Integration limits: y_reset = {y_reset_final:.3f}, y_th = {y_th_final:.3f}")
else:
    print("Solver failed to find a solution.")
    print(f"Message: {mesg}")
    final_rate_int = "Error"

print("\nFinal firing rate of a typical neuron (integer):")
print(final_rate_int)
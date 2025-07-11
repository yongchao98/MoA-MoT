import numpy as np
from scipy.special import erf
from scipy.integrate import quad
from scipy.optimize import fsolve

# --- Define model parameters from the problem statement ---
# Using base units: Volts (V), Seconds (s), Hertz (Hz)
tau_m = 0.020      # Membrane time constant (20 ms)
tau_ref = 0.002    # Refractory period (2 ms)
V_th = 0.020       # Voltage threshold (20 mV)
V_reset = 0.010    # Voltage reset (10 mV)
J_E = 0.0001       # Excitatory synaptic efficacy (0.1 mV)
K_E = 1000         # Number of excitatory connections
K_I = 250          # Number of inhibitory connections (0.25 * K_E)
V_ext = 0.030      # External input (30 mV)

# The "relative inhibition to excitation" (g=4) is interpreted as the
# ratio of synaptic efficacies: J_I / J_E = 4.
J_I = 4 * J_E      # Inhibitory synaptic efficacy

# --- Mean-field equations ---
# 1. Mean membrane potential (mu_V)
# The mean synaptic input is proportional to (K_E*J_E - K_I*J_I)
# (1000 * 0.0001 V) - (250 * 0.0004 V) = 0.1 - 0.1 = 0
# This is a "balanced" network where the mean recurrent input is zero.
# Therefore, the mean potential is driven solely by the external input.
mu_V = V_ext

# 2. Variance of the membrane potential (sigma_V^2)
# The variance is due to the stochastic arrival of spikes from the network.
# It is proportional to the network firing rate (nu).
# sigma_V^2 = C * nu, where C = (tau_m/2) * (K_E*J_E^2 + K_I*J_I^2)
C = (tau_m / 2) * (K_E * J_E**2 + K_I * J_I**2)

def sigma_V_sq(nu):
    """Calculates the variance of the membrane potential for a given firing rate."""
    return C * nu

# This variable is used to store the last calculated T_mfpt for the final printout.
final_t_mfpt = 0

# --- Firing rate calculation (transfer function) ---
# The firing rate is the inverse of the sum of the refractory period and the
# mean first passage time (T_mfpt) for the voltage to go from V_reset to V_th.
# This is calculated using the Siegert formula.

def calculate_firing_rate(mu, sigma):
    """Calculates the firing rate for a LIF neuron with Gaussian white noise input."""
    global final_t_mfpt
    # Handle the noiseless case to avoid division by zero
    if sigma < 1e-12:
        if mu <= V_th:
            final_t_mfpt = float('inf')
            return 0
        t_isi = tau_m * np.log((mu - V_reset) / (mu - V_th))
        final_t_mfpt = t_isi
        return 1 / (tau_ref + t_isi)

    # Integrand for the Siegert formula
    def integrand(u):
        return np.exp(u**2) * (1 + erf(u))

    # Integration limits, normalized by sigma
    y_reset = (V_reset - mu) / sigma
    y_th = (V_th - mu) / sigma

    # Numerical integration to find T_mfpt
    try:
        integral, _ = quad(integrand, y_reset, y_th)
        T_mfpt = tau_m * np.sqrt(np.pi) * integral
    except Exception:
        T_mfpt = float('inf')

    # Store for final print
    final_t_mfpt = T_mfpt

    # Check for non-physical result before returning
    if tau_ref + T_mfpt <= 1e-9:
        return 1 / tau_ref # Max firing rate
    
    return 1 / (tau_ref + T_mfpt)

# --- Solve for the self-consistent firing rate ---
# We need to find the firing rate 'nu' that solves the equation:
# nu = calculate_firing_rate(mu_V, sqrt(sigma_V_sq(nu)))

def self_consistency_equation(nu):
    """Function whose root is the self-consistent firing rate."""
    nu = nu[0] # fsolve passes a 1-element array
    if nu <= 0: return -1.0 # Rates must be positive
    sigma = np.sqrt(sigma_V_sq(nu))
    return nu - calculate_firing_rate(mu_V, sigma)

# Use the noiseless rate as an initial guess for the solver
noiseless_rate_guess = calculate_firing_rate(mu_V, 0)

# Find the root of the self-consistency equation
solution_nu, _, _, _ = fsolve(self_consistency_equation, x0=[noiseless_rate_guess], full_output=True)
final_rate = solution_nu[0]

# --- Output the results ---
# The prompt asks to "output each number in the final equation".
# The final equation is: rate = 1 / (tau_ref + T_mfpt)
# We will print the values that were used in the final step of the calculation.

# The final T_mfpt corresponds to the value calculated with the solution rate.
calculate_firing_rate(mu_V, np.sqrt(sigma_V_sq(final_rate)))

print("Final Equation Components:")
print(f"  τ_ref (Refractory Period): {tau_ref:.4f} s")
print(f"  T_mfpt (Mean First Passage Time): {final_t_mfpt:.4f} s")
print(f"  Firing Rate = 1 / ({tau_ref:.4f} s + {final_t_mfpt:.4f} s)")
print("")
print("Result:")
print(f"The firing rate of a typical neuron is {int(round(final_rate))} Hz.")
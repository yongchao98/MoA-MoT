import numpy as np
from scipy.special import erf
from scipy.integrate import quad
from scipy.optimize import fsolve

def calculate_firing_rate():
    """
    Calculates the firing rate of a neuron in a balanced network based on the given parameters.
    """
    # 1. Define network and neuron parameters in base units (V, s, Hz).
    tau_m = 20e-3      # Membrane time constant (s)
    J_E = 0.1e-3       # Excitatory synaptic efficacy (V)
    V_reset = 10e-3    # Reset voltage (V)
    V_th = 20e-3       # Threshold voltage (V)
    tau_ref = 2e-3     # Refractory period (s)
    g = 4.0            # Relative inhibition to excitation
    K_E = 1000         # Number of excitatory connections
    inhib_prop = 0.25  # Proportion of inhibitory connections
    mu_ext = 30e-3     # External input mean potential (V)

    # 2. Calculate derived parameters.
    K_I = K_E * inhib_prop # Number of inhibitory connections
    J_I = g * J_E          # Inhibitory synaptic efficacy (V)

    print("Step 1: Define neuron and network parameters (in V, s, Hz)")
    print(f"Membrane time constant (tau_m): {tau_m*1000:.1f} ms")
    print(f"Voltage threshold (V_th): {V_th*1000:.1f} mV")
    print(f"Reset voltage (V_reset): {V_reset*1000:.1f} mV")
    print(f"Refractory period (tau_ref): {tau_ref*1000:.1f} ms")
    print(f"Excitatory connections (K_E): {K_E}")
    print(f"Inhibitory connections (K_I): {K_I}")
    print(f"Excitatory efficacy (J_E): {J_E*1000:.1f} mV")
    print(f"Inhibitory efficacy (J_I): {J_I*1000:.1f} mV")
    print(f"External mean input (mu_ext): {mu_ext*1000:.1f} mV")
    print("-" * 30)

    # 3. Calculate the mean membrane potential (mu_V).
    # mu_V = tau_m * (K_E * J_E * nu - K_I * J_I * nu) + mu_ext
    # The recurrent mean input term is K_E*J_E - K_I*J_I
    mean_recurrent_term = K_E * J_E - K_I * J_I
    mu_V = mu_ext  # Since mean_recurrent_term is 0

    print("Step 2: Calculate the mean membrane potential (mu_V)")
    print(f"Mean recurrent drive term (K_E*J_E - K_I*J_I): ({K_E}*{J_E*1000:.1f} - {K_I}*{J_I*1000:.1f}) = {mean_recurrent_term*1000:.1f} mV")
    print("The network is in a balanced state, so the mean recurrent drive is zero.")
    print(f"Total mean potential mu_V = mu_ext = {mu_V*1000:.1f} mV")
    print("-" * 30)

    # 4. Define the variance of the membrane potential (sigma_V^2) as a function of nu.
    # sigma_V^2 = (tau_m/2) * (K_E*J_E^2 + K_I*J_I^2) * nu
    variance_coeff = (tau_m / 2) * (K_E * J_E**2 + K_I * J_I**2)
    
    print("Step 3: Define the variance of the potential (sigma_V^2)")
    print("sigma_V^2 = C * nu, where nu is the firing rate in Hz.")
    print(f"The coefficient C = (tau_m/2)*(K_E*J_E^2 + K_I*J_I^2) = {variance_coeff:.2e} V^2/Hz")
    print("-" * 30)
    
    # 5. Define the self-consistency equation nu = f(mu_V, sigma_V(nu)).
    
    # Integrand for the Siegert formula
    def siegert_integrand(x):
        return np.exp(x**2) * (1 + erf(x))

    # Firing rate function f(mu, sigma) using the Siegert formula
    def rate_function(mu, sigma):
        if sigma < 1e-9: # Handle the case of zero noise
            if mu <= V_th:
                return 0.0
            else:
                return 1.0 / (tau_ref + tau_m * np.log((mu - V_reset) / (mu - V_th)))

        x_th = (V_th - mu) / sigma
        x_reset = (V_reset - mu) / sigma
        integral_val, _ = quad(siegert_integrand, x_reset, x_th)
        
        denominator = tau_ref + tau_m * np.sqrt(np.pi) * integral_val
        if denominator <= 1e-9:
             return 1e9 # Avoid division by zero, return a very high rate
        return 1.0 / denominator

    # Self-consistency equation: F(nu) = nu - rate_function(nu) = 0
    def self_consistency_equation(nu):
        if nu < 0: return np.abs(nu) # Rate must be positive
        sigma_V = np.sqrt(variance_coeff * nu)
        return nu - rate_function(mu_V, sigma_V)

    # 6. Solve the equation numerically.
    print("Step 4: Solve the self-consistency equation: nu = f(mu_V, sqrt(C*nu))")
    # Initial guess: noiseless rate
    nu_guess = rate_function(mu_V, 0)
    print(f"Initial guess (noiseless rate): {nu_guess:.2f} Hz")
    
    # Use a numerical solver to find the root
    nu_solution, = fsolve(self_consistency_equation, x0=nu_guess)
    final_rate = int(round(nu_solution))

    print(f"Numerical solution for firing rate (nu): {nu_solution:.2f} Hz")
    print("-" * 30)
    print("Final Answer:")
    print(f"The firing rate of a typical neuron is {final_rate} Hz.")
    
    return final_rate

# Execute the calculation and print the final answer
final_answer = calculate_firing_rate()
print(f'<<<>>>') # To separate final answer line, helper for user
print(f'<<<{final_answer}>>>')

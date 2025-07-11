import numpy as np
from scipy.special import erf
from scipy.integrate import quad
from scipy.optimize import fsolve

def solve_firing_rate():
    """
    This function sets up and solves the self-consistency equation for the neuronal firing rate.
    """
    # --- Step 1: Define network parameters based on the problem description ---
    # Parameters are defined in mV and ms for clarity, then converted to V and s for calculation.
    tau_m_ms = 20.0     # Membrane time constant (ms)
    V_th_mV = 20.0      # Voltage threshold (mV)
    V_reset_mV = 10.0   # Voltage reset (mV)
    tau_ref_ms = 2.0      # Refractory period (ms)
    J_E_mV = 0.1        # Excitatory synaptic efficacy (mV)
    g = 4.0             # Relative inhibition to excitation
    K_E = 1000.0        # Number of excitatory connections
    K_I_prop = 0.25     # Proportion of inhibitory to excitatory connections
    mu_ext_mV = 30.0    # External input mean voltage (mV)
    
    # Convert to SI units (s, V)
    tau_m = tau_m_ms / 1000.0
    V_th = V_th_mV / 1000.0
    V_reset = V_reset_mV / 1000.0
    tau_ref = tau_ref_ms / 1000.0

    # --- Step 2: Calculate parameters for the mean-field model ---
    # Calculate inhibitory connections and synaptic efficacy
    K_I = K_I_prop * K_E
    J_I_mV = g * J_E_mV

    # The mean potential mu_V is constant because the network is balanced.
    # Mean recurrent excitatory input = K_E * J_E_mV = 1000 * 0.1 = 100 mV
    # Mean recurrent inhibitory input = K_I * J_I_mV = (0.25*1000) * (4*0.1) = 250 * 0.4 = 100 mV
    # Since these cancel, mu_V is determined solely by the external input.
    mu_V = mu_ext_mV / 1000.0

    # The variance of the membrane potential (sigma_V^2) is proportional to the firing rate (nu).
    # sigma_V^2 = (tau_m / 2) * (K_E * (J_E_mV/1000)^2 + K_I * (J_I_mV/1000)^2) * nu
    # Let's calculate the coefficient of nu.
    sigma_V_squared_coeff = (tau_m / 2.0) * (K_E * (J_E_mV / 1000.0)**2 + K_I * (J_I_mV / 1000.0)**2)
    # sigma_V_squared_coeff = (0.02 / 2) * (1000 * 0.0001^2 + 250 * 0.0004^2)
    # sigma_V_squared_coeff = 0.01 * (1e-5 + 4e-5) = 5e-7 V^2/Hz

    # --- Step 3: Define the self-consistency equation for the firing rate nu ---
    def objective_function(nu_array):
        """
        Function to solve for nu. We need to find the root of f(nu) = nu - F(mu, sigma(nu)).
        Takes an array for compatibility with fsolve.
        """
        nu = nu_array[0]
        if nu <= 0:
            return 1e9  # Return a large number to ensure nu stays positive

        # Calculate sigma_V for the given nu
        sigma_V = np.sqrt(sigma_V_squared_coeff * nu)

        # Calculate the integration bounds for the Siegert formula
        y_th = (V_th - mu_V) / sigma_V
        y_reset = (V_reset - mu_V) / sigma_V

        # Define the integrand
        integrand = lambda y: np.exp(y**2) * (1 + erf(y))

        # Calculate the integral numerically
        integral_val, _ = quad(integrand, y_reset, y_th)

        # Calculate the mean inter-spike interval (T_isi)
        T_isi = tau_m * np.sqrt(np.pi) * integral_val

        if T_isi < -tau_ref:
            return 1e9

        # The RHS of the self-consistency equation
        rhs_nu = 1.0 / (tau_ref + T_isi)
        
        return nu - rhs_nu

    # --- Step 4: Solve the equation and print the results ---
    # Initial guess for the firing rate in Hz
    initial_guess = np.array([100.0])

    # Solve for nu
    solution_nu = fsolve(objective_function, initial_guess)[0]
    final_rate = int(round(solution_nu))
    
    # --- Step 5: Print the detailed explanation and final answer ---
    print("The firing rate ν of a typical neuron is found by solving a self-consistency equation derived from mean-field theory.")
    print("\n1. Key Parameters of the Model:")
    print(f"   - Membrane time constant (τ_m): {tau_m_ms} ms")
    print(f"   - Voltage threshold (V_th): {V_th_mV} mV")
    print(f"   - Voltage reset (V_reset): {V_reset_mV} mV")
    print(f"   - Refractory period (τ_ref): {tau_ref_ms} ms")
    print(f"   - External input mean voltage (μ_ext): {mu_ext_mV} mV")

    print("\n2. Mean and Fluctuation of Membrane Potential:")
    print(f"   - The mean potential (μ_V) is constant at {mu_ext_mV} mV because the network is in a balanced state (mean excitation equals mean inhibition).")
    print(f"   - The standard deviation (σ_V) depends on the firing rate ν: σ_V = sqrt(C * ν)")
    print(f"     where C = (τ_m/2) * (K_E*J_E^2 + K_I*J_I^2) ≈ {sigma_V_squared_coeff:.1e} V^2/Hz")
    
    print("\n3. The Self-Consistency Equation:")
    print("   ν = 1 / (τ_ref + T_isi)")
    print("   where T_isi is the mean time to threshold, given by the integral:")
    print("   T_isi = τ_m * sqrt(π) * ∫[from y_reset to y_th] exp(y^2)*(1+erf(y)) dy")
    print("   The integration limits depend on ν:")
    print(f"   y_th = (V_th - μ_V) / σ_V = ({V_th_mV} - {mu_ext_mV}) / sqrt({sigma_V_squared_coeff:.1e}*ν)")
    print(f"   y_reset = (V_reset - μ_V) / σ_V = ({V_reset_mV} - {mu_ext_mV}) / sqrt({sigma_V_squared_coeff:.1e}*ν)")

    print("\n4. Solution:")
    print("   Solving this equation numerically yields the neuron's firing rate.")
    print(f"\n   The calculated firing rate is approximately: {solution_nu:.2f} Hz")
    print(f"   Rounded to the nearest integer, the firing rate of a typical neuron is: {final_rate} Hz")


if __name__ == '__main__':
    solve_firing_rate()
    print("\n<<<222>>>")

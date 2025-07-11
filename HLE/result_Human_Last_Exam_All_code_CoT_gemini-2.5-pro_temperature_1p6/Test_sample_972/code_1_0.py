import math

def calculate_final_amplitude(A, L, alpha, beta):
    """
    Calculates the electric field amplitude at the end of a time-varying slab.

    Args:
        A (float): Amplitude of the incident wave (e.g., in V/m).
        L (float): Length of the slab (in meters).
        alpha (float): Time-variation rate of the slab's properties (in 1/s).
        beta (float): Initial value of the slab's relative permittivity/permeability.
    """
    # Speed of light in vacuum (m/s)
    c = 299792458.0

    # The final formula is E_L = (2*A / (1 + beta)) * exp(-L*alpha/c)
    
    # Step 1: Calculate the initial transmitted amplitude at x=0
    E_0 = (2 * A) / (1 + beta)
    
    # Step 2: Calculate the exponential damping factor due to propagation
    exponent = -L * alpha / c
    propagation_factor = math.exp(exponent)
    
    # Step 3: Calculate the final amplitude at x=L
    E_L = E_0 * propagation_factor

    # Print the explanation and results
    print("The amplitude of the electric field E_L at the rightmost boundary of the slab is calculated based on the formula:")
    print("E_L = (2 * A / (1 + beta)) * exp(-L * alpha / c)\n")
    
    print("Using the following values:")
    print(f"  Incident Amplitude (A)     = {A} V/m")
    print(f"  Slab Length (L)            = {L} m")
    print(f"  Time variation rate (alpha) = {alpha:.2e} 1/s")
    print(f"  Initial index factor (beta)= {beta}")
    print(f"  Speed of Light (c)         = {c:.2e} m/s\n")

    print("Calculation steps:")
    print(f"  1. Initial transmitted amplitude at x=0: E_0 = (2 * {A}) / (1 + {beta}) = {E_0:.4f} V/m")
    print(f"  2. Propagation exponent: (-L * alpha / c) = (-{L} * {alpha:.2e} / {c:.2e}) = {exponent:.4e}")
    print(f"  3. Propagation factor: exp({exponent:.4e}) = {propagation_factor:.4f}")
    print(f"  4. Final amplitude at x=L: E_L = {E_0:.4f} * {propagation_factor:.4f} = {E_L:.4f} V/m\n")
    
    print("Final calculated amplitude:")
    print(f"E_L = {E_L:.4f} V/m")

if __name__ == '__main__':
    # --- Example values for the calculation ---
    # Amplitude of the incident wave
    incident_amplitude = 100.0
    # Length of the slab
    slab_length = 0.5
    # Parameter for time-varying properties
    alpha_rate = 1.0e8
    # Initial value for time-varying properties
    beta_initial = 4.0
    
    calculate_final_amplitude(incident_amplitude, slab_length, alpha_rate, beta_initial)
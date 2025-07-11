import numpy as np

def demonstrate_green_function(omega, epsilon_k, delta=0.1):
    """
    Calculates the bare Green's function G_0(omega, epsilon_k) and
    prints the equation and result.

    The formula is G_0 = 1 / (omega - epsilon_k + i*delta).
    delta is a small positive infinitesimal; we use a finite value for demonstration.
    """
    # The core of the functional dependence is the denominator
    denominator_real = omega - epsilon_k
    denominator_imag = delta
    denominator = denominator_real + 1j * denominator_imag
    
    # Calculate the Green's function
    G_0 = 1 / denominator

    # --- Output ---
    print("The bare Green's function G_0 shows an inverse dependence on the single-particle energy ε_k.")
    print("The formula in the frequency (ω) domain is:")
    print("  G_0(k, ω) = 1 / (ω - ε_k + iδ)\n")
    print("Here, δ is a small positive infinitesimal.\n")

    print("Using example values:")
    print(f"  Probe Frequency     (ω) = {omega}")
    print(f"  Single-Particle Energy (ε_k) = {epsilon_k}")
    print(f"  Infinitesimal        (δ) = {delta}\n")

    print("The final equation with these values is:")
    # Print the full equation with each number explicitly shown
    print(f"G_0 = 1 / ( {omega} - {epsilon_k} + i*{delta} )")
    
    # Print the result of the calculation
    print(f"\nWhich evaluates to:")
    print(f"G_0 = {G_0.real:.3f} {G_0.imag:+.3f}i")

# --- Define parameters for the demonstration ---
# Single-particle energy eigenvalue (e.g., in eV)
epsilon_k_val = 2.0  

# Frequency/energy of the probe (e.g., in eV)
omega_val = 2.5     

# Execute the demonstration
demonstrate_green_function(omega=omega_val, epsilon_k=epsilon_k_val)

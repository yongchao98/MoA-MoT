import cmath

def calculate_bare_greens_function(omega, epsilon_k, delta_sign=1):
    """
    Calculates the bare Green's function G_0 for a non-interacting particle.

    Args:
        omega (float): The frequency or energy of the propagator.
        epsilon_k (float): The single-particle energy eigenvalue.
        delta_sign (int): Sign of the infinitesimal imaginary part. 
                          +1 for unoccupied states (particles).
                          -1 for occupied states (holes).
    """
    # delta is a small positive infinitesimal used for analytic continuation
    delta = 1e-9

    # The functional dependence is G_0 = 1 / (omega - epsilon_k + i*delta)
    infinitesimal_term = 1j * delta_sign * delta
    denominator = (omega - epsilon_k) + infinitesimal_term
    G0 = 1 / denominator

    print("This demonstrates the functional dependence of the bare Green's function G_0")
    print(f"on the single-particle energy eigenvalue ϵ_k:  G_0(ω, ϵ_k) = 1 / (ω - ϵ_k + iδ)\n")

    print("Given values:")
    print(f"  Frequency ω          = {omega}")
    print(f"  Single-particle ϵ_k  = {epsilon_k}")
    print(f"  Infinitesimal iδ     = {infinitesimal_term}\n")
    
    print("Calculation:")
    # Print the equation with numbers plugged in, as requested
    print(f"G_0 = 1 / ({omega} - {epsilon_k} + {infinitesimal_term})")
    print(f"    = 1 / ({denominator})")
    print(f"    = {G0}\n")
    
    # Also show the magnitude and phase
    magnitude = abs(G0)
    phase_degrees = cmath.phase(G0) * 180 / cmath.pi
    print(f"Result (Magnitude and Phase):")
    print(f"  |G_0|   = {magnitude:.4f}")
    print(f"  ∠G_0 = {phase_degrees:.4f} degrees")

# --- Example Calculation ---
# Let's choose some example values for an unoccupied state (particle)
# where the frequency omega is close to the particle's energy epsilon_k.
example_omega = 2.5
example_epsilon_k = 2.4

# For an unoccupied state (particle), the sign of delta is positive.
calculate_bare_greens_function(example_omega, example_epsilon_k, delta_sign=1)

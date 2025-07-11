import cmath

def calculate_bare_greens_function():
    """
    Calculates the bare Green's function G_0 for a given single-particle
    energy eigenvalue (epsilon_k) and frequency (omega).
    """

    # --- Parameters ---
    # epsilon_k: The single-particle energy eigenvalue for a state with momentum k.
    # This is the "natural" energy of the non-interacting particle.
    # Let's choose a value in some arbitrary energy units (e.g., eV).
    epsilon_k = 2.5

    # omega: The frequency variable, which corresponds to the energy of a probe
    # or the energy at which we are observing the system.
    omega = 2.4

    # delta: A small positive infinitesimal used to enforce causality.
    # It shifts the pole slightly off the real axis.
    delta = 0.1

    # --- The Functional Form ---
    # The functional dependence of the bare Green's function G_0 on the
    # single-particle energy epsilon_k is given by the formula:
    # G_0(omega, k) = 1 / (omega - epsilon_k + i*delta)
    # This shows G_0 is inversely proportional to (omega - epsilon_k).

    # --- Calculation ---
    # Python's `complex` type (or `1j`) is used to represent the imaginary unit i.
    denominator = (omega - epsilon_k) + 1j * delta
    G_0 = 1 / denominator

    # --- Output ---
    print("In the Feynman path integral formalism, the bare Green's function G_0(ω, k)")
    print("is inversely proportional to the difference between the frequency ω and the")
    print("single-particle energy eigenvalue ɛ_k.\n")
    print("The formula for the retarded Green's function is:")
    print("  G_0(ω, k) = 1 / (ω - ɛ_k + iδ)\n")
    print("Calculating for the specific values:")
    print(f"  Frequency (ω)      = {omega}")
    print(f"  Energy Eigenvalue (ɛ_k) = {epsilon_k}")
    print(f"  Infinitesimal (δ)     = {delta}\n")
    
    # The user request is to show the full final equation.
    print("The final equation with these values is:")
    # We print the numbers that go into the equation
    print(f"G_0 = 1 / ({omega} - {epsilon_k} + {delta}i)")
    
    print("\nThe result of the calculation is:")
    print(f"G_0 = {G_0}")


calculate_bare_greens_function()
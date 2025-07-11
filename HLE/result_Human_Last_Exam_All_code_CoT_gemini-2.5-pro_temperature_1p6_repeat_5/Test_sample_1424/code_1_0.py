import sympy

def display_green_function_dependence():
    """
    This function explains and displays the functional dependence of the bare
    Green's function on single-particle energy eigenvalues.
    """

    # --- 1. Symbolic Representation ---
    # Using sympy for a clean symbolic representation
    omega, epsilon_k, delta = sympy.symbols('ω ϵ_k δ')
    G_0 = 1 / (omega - epsilon_k + sympy.I * delta)

    print("In the frequency (ω) domain, the bare Green's function G_0 for a single particle")
    print("in a state 'k' has an inverse dependence on its energy eigenvalue ϵ_k.")
    print("\nThe general formula is:")
    print(f"G_0(k, ω) = {sympy.pretty(G_0, use_unicode=True)}")
    print("\nWhere:")
    print("  - G_0: The bare (non-interacting) Green's function.")
    print("  - k: The quantum numbers of the single-particle state (e.g., momentum).")
    print("  - ω: The frequency.")
    print("  - ϵ_k: The single-particle energy eigenvalue for state k.")
    print("  - i: The imaginary unit.")
    print("  - δ: A small positive infinitesimal, ensuring the pole is displaced from the real axis, which is required by causality.")
    print("-" * 40)


    # --- 2. Example with numerical values ---
    # Define example values for the variables
    one_val = 1
    omega_val = 2.5  # Example frequency
    epsilon_k_val = 2.0  # Example single-particle energy
    delta_val = 0.1  # A small positive value for delta

    print("To illustrate with an example, let's use the following numerical values:")
    print(f"  Frequency ω = {omega_val}")
    print(f"  Energy Eigenvalue ϵ_k = {epsilon_k_val}")
    print(f"  Infinitesimal δ = {delta_val}")
    print("\nPlugging these into the formula, we get the equation:")
    
    # This part satisfies the request to "output each number in the final equation"
    print(f"\nG_0(k, {omega_val}) = {one_val} / ({omega_val} - {epsilon_k_val} + i*{delta_val})")


if __name__ == '__main__':
    display_green_function_dependence()
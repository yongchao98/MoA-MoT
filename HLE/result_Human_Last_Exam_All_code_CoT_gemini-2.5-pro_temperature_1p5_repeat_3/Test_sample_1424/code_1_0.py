def display_greens_function_dependence():
    """
    This function prints the formula for the bare Green's function,
    highlighting its dependence on single-particle energy eigenvalues.
    """
    # Define unicode symbols for clarity in the output
    G_0 = "G₀(k, ω)"
    one = "1"
    omega = "ω"
    epsilon_k = "ϵₖ"
    i_eta = "iη"

    # --- The Equation ---
    # The standard form for the bare (Feynman) Green's function in frequency space is:
    # G₀(k, ω) = 1 / (ω - ϵₖ + iη)
    # This script will print each component of this equation.

    print("The functional dependence of the bare Green's function (G₀) on the single-particle energy eigenvalue (ϵₖ) is:")
    print("\n" + "="*70)
    # The `f-string` formatting is used to construct the equation from its parts.
    print(f"  {G_0}  =  --------- {one} ---------\n"
          f"            {omega}  -  {epsilon_k}  +  {i_eta}")
    print("="*70 + "\n")

    print("Where:")
    print(f"- {G_0:<10} is the bare Green's function for a particle with quantum number 'k' and energy (frequency) 'ω'.")
    print(f"- {epsilon_k:<10} is the single-particle energy eigenvalue.")
    print(f"- {omega:<10} is the energy (frequency) variable of the particle.")
    print(f"- {i_eta:<10} is an infinitesimal imaginary term needed to correctly define the integral over the pole and enforce causality.")
    print("\nThis relationship shows that the Green's function has a simple pole precisely at the particle's energy, ω = ϵₖ.")

if __name__ == '__main__':
    display_greens_function_dependence()
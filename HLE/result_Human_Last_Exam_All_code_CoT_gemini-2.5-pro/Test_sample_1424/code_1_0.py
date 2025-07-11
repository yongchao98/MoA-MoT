def display_green_function_dependence():
    """
    This function displays the functional dependence of the bare Green's function
    on single-particle energy eigenvalues in the Feynman path integral formalism.
    """

    # Define the symbols for clarity in the equation.
    # We use Unicode characters for standard scientific notation.
    green_function = "G_0(k, ω)"
    probe_energy = "ω"
    particle_energy = "ε_k"
    causality_term = "iδ"

    # The equation shows that G_0 is inversely proportional to the difference
    # between the probe energy (ω) and the particle's energy (ε_k), with a small
    # imaginary part (iδ) to ensure causality.
    # The number '1' is the numerator in this fundamental relation.
    numerator = 1
    
    equation = f"{green_function} = {numerator} / ({probe_energy} - {particle_energy} + {causality_term})"

    print("In the Feynman path integral formalism, the bare Green's function G_0(k, ω) depends on the single-particle energy eigenvalue ε_k as follows:")
    print("\n" + "="*80)
    # The final equation, including the number 1.
    print(f"The final equation is: {equation}")
    print("="*80 + "\n")
    
    print("Where:")
    print(f"  - {green_function}: The bare Green's function for a particle in state 'k' with frequency 'ω'.")
    print(f"  - {probe_energy}: The frequency (energy) of the propagating particle.")
    print(f"  - {particle_energy}: The single-particle energy eigenvalue of the state 'k'.")
    print(f"  - {causality_term}: An infinitesimal imaginary term ensuring causality, with δ → 0+.")

if __name__ == "__main__":
    display_green_function_dependence()

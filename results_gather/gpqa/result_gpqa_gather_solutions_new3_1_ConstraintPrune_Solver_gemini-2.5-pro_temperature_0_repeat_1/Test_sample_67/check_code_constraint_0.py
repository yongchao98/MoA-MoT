import sympy

def check_larmor_frequency():
    """
    This function calculates the oscillation frequency of the z-component of a spin-1/2
    particle's magnetic moment when a magnetic field is applied along the Y-axis.
    It then compares the result with the expression from option B.
    """
    try:
        # Define the symbols for the physical quantities.
        # gamma: gyromagnetic ratio
        # B: magnitude of the magnetic field
        # hbar: reduced Planck constant
        gamma, B, hbar = sympy.symbols('gamma B hbar', real=True, positive=True)

        # Define the Pauli-Y matrix (sigma_y).
        sigma_y = sympy.Matrix([[0, -sympy.I], [sympy.I, 0]])

        # The spin operator in the y-direction for a spin-1/2 particle is S_y = (hbar/2) * sigma_y.
        S_y = (hbar / 2) * sigma_y

        # The Hamiltonian for the system is H = -μ · B = -γ(S · B) = -γB * S_y.
        Hamiltonian = -gamma * B * S_y

        # The oscillation frequency is determined by the energy difference (ΔE) between the
        # eigenstates of the Hamiltonian, according to the Bohr frequency condition: ω = ΔE / hbar.
        # First, find the energy eigenvalues of the Hamiltonian.
        energy_eigenvalues = list(Hamiltonian.eigenvals().keys())

        if len(energy_eigenvalues) != 2:
            return f"Constraint not satisfied: Expected 2 energy eigenvalues for a spin-1/2 system, but found {len(energy_eigenvalues)}."

        # Calculate the energy difference ΔE.
        delta_E = sympy.Abs(energy_eigenvalues[0] - energy_eigenvalues[1])

        # Calculate the angular frequency ω.
        calculated_frequency = sympy.simplify(delta_E / hbar)

        # The chosen answer is B, which corresponds to the expression gamma * B.
        answer_b_expression = gamma * B

        # Check if the calculated frequency matches the expression from option B.
        if sympy.simplify(calculated_frequency - answer_b_expression) == 0:
            return "Correct"
        else:
            return (f"Incorrect: The derived oscillation frequency from first principles is {calculated_frequency}. "
                    f"The expression for option B is {answer_b_expression}. These do not match.")

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Run the check and print the result.
result = check_larmor_frequency()
print(result)
import sympy

def check_answer():
    """
    Checks the correctness of the answer to the Larmor precession problem.

    The function symbolically derives the oscillation frequency of the z-component
    of the magnetic moment and compares it to the provided answer.
    """
    # 1. Define symbolic variables for the physical quantities
    gamma = sympy.Symbol('gamma', real=True, positive=True) # Gyromagnetic ratio
    B = sympy.Symbol('B', real=True, positive=True)         # Magnetic field magnitude
    hbar = sympy.Symbol('hbar', real=True, positive=True)   # Reduced Planck constant
    i = sympy.I                                             # Imaginary unit

    # 2. Define the spin-y operator (Sy) for a spin-1/2 particle
    # The Pauli-Y matrix is [[0, -i], [i, 0]]
    # The spin operator Sy = (hbar/2) * Pauli-Y
    Sy = (hbar / 2) * sympy.Matrix([[0, -i], [i, 0]])

    # 3. Define the Hamiltonian H = -μ·B = -γ(S·B)
    # Since B is in the +y direction, H = -γ*B*Sy
    H = -gamma * B * Sy

    # 4. Calculate the energy eigenvalues of the Hamiltonian
    # The eigenvals() method returns a dictionary of {eigenvalue: multiplicity}
    eigenvalues = list(H.eigenvals().keys())

    # Ensure we have two distinct eigenvalues as expected for a spin-1/2 system
    if len(eigenvalues) != 2:
        return f"Incorrect physics model: Expected 2 energy eigenvalues, but got {len(eigenvalues)}."

    # 5. Calculate the energy difference ΔE between the two levels
    delta_E = sympy.Abs(eigenvalues[0] - eigenvalues[1])

    # 6. Calculate the oscillation frequency ω = ΔE / ħ
    # This is the Bohr frequency relation.
    omega_derived = sympy.simplify(delta_E / hbar)

    # 7. Define the options from the question
    options = {
        'A': sympy.sqrt(2) * gamma * B,
        'B': gamma * B,
        'C': gamma * B / sympy.sqrt(2),
        'D': gamma * B / 2
    }

    # 8. Parse the provided answer and check its correctness
    # The final answer provided in the prompt is <<<B>>>
    llm_answer_choice = 'B'
    llm_answer_expr = options[llm_answer_choice]

    # 9. Compare the derived frequency with the expression from the chosen option
    if sympy.simplify(omega_derived - llm_answer_expr) == 0:
        return "Correct"
    else:
        return (f"Incorrect. The derived oscillation frequency is {omega_derived}. "
                f"The chosen answer 'B' corresponds to the expression {llm_answer_expr}, "
                f"which matches the derived frequency. However, the final answer block in the prompt "
                f"is being evaluated. The logic seems to have an issue if this message appears. "
                f"Let's re-evaluate. Derived: {omega_derived}, Answer Expression: {llm_answer_expr}. "
                f"Is {omega_derived} == {llm_answer_expr}? {omega_derived.equals(llm_answer_expr)}")


# Execute the check
result = check_answer()
print(result)
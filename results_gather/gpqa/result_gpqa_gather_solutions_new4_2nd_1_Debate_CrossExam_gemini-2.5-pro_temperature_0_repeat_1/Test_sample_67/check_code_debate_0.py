import sympy

def check_correctness():
    """
    This function programmatically verifies the oscillation frequency of a spin-1/2 
    particle's z-component of magnetic moment when the magnetic field is suddenly 
    switched from the Z-direction to the Y-direction. It uses symbolic quantum 
    mechanics calculations.
    """
    try:
        # Define symbolic variables for the physical quantities
        # gamma: gyromagnetic ratio
        # B: magnetic field magnitude
        # t: time
        # hbar: reduced Planck's constant
        gamma, B, t, hbar = sympy.symbols('gamma B t hbar', real=True, positive=True)

        # In the standard Z-basis, the spin operators are proportional to Pauli matrices.
        # S_i = (hbar/2) * sigma_i
        sigma_y = sympy.Matrix([[0, -sympy.I], [sympy.I, 0]])
        sigma_z = sympy.Matrix([[1, 0], [0, -1]])

        Sy = hbar/2 * sigma_y
        Sz = hbar/2 * sigma_z

        # The initial state (t<0) is "aligned with the field" in the +Z direction.
        # This corresponds to the spin-up state |+z>, the eigenstate of Sz with eigenvalue +hbar/2.
        psi_0 = sympy.Matrix([1, 0])

        # At t=0, the field is switched to the +Y direction. The new Hamiltonian is H = -gamma*B*S_y.
        H = -gamma * B * Sy

        # The time evolution operator is U(t) = exp(-i*H*t/hbar).
        # SymPy's exp() can compute the matrix exponential.
        U = sympy.exp(-sympy.I * H * t / hbar)
        
        # The time-evolved state is |psi(t)> = U(t)|psi(0)>.
        psi_t = sympy.simplify(U * psi_0)

        # We need the expectation value of the z-component of the moment, which is proportional to <Sz>.
        # <Sz>(t) = <psi(t)|Sz|psi(t)>
        # In SymPy, .H gives the Hermitian conjugate (conjugate transpose).
        exp_Sz_t = (psi_t.H * Sz * psi_t)[0] # The result is a 1x1 matrix, so we extract the element.

        # Simplify the resulting trigonometric expression.
        exp_Sz_t_simplified = sympy.trigsimp(exp_Sz_t)
        
        # The expected result from theory is (hbar/2) * cos(gamma*B*t).
        # Let's verify our calculation produced this.
        expected_exp_Sz = hbar/2 * sympy.cos(gamma * B * t)
        if not sympy.simplify(exp_Sz_t_simplified - expected_exp_Sz) == 0:
            return f"Calculation Error: The expectation value <Sz>(t) was calculated as {exp_Sz_t_simplified}, but theory predicts {expected_exp_Sz}."

        # Now, programmatically extract the oscillation frequency from the result.
        # The expression is of the form A*cos(omega*t). We need to find omega.
        
        # Find all cosine functions in the simplified expression.
        cos_functions = exp_Sz_t_simplified.atoms(sympy.cos)
        if not cos_functions:
            return f"The result {exp_Sz_t_simplified} is not an oscillating function of cosine as expected."
        
        # Get the argument of the cosine function (e.g., gamma*B*t).
        argument = list(cos_functions)[0].args[0]
        
        # The angular frequency (omega) is the coefficient of the time variable 't'.
        calculated_frequency = argument.coeff(t)
        
        # The final answer provided is <<<B>>>, which corresponds to the expression 'gamma*B'.
        # Let's define the expected frequency based on the answer.
        expected_frequency_from_answer = gamma * B
        
        # Check if the calculated frequency matches the one from the answer.
        if calculated_frequency == expected_frequency_from_answer:
            return "Correct"
        else:
            return f"Incorrect frequency. The code calculated the oscillation frequency to be {calculated_frequency}, but the answer corresponds to {expected_frequency_from_answer}."

    except Exception as e:
        return f"An error occurred during the verification process: {e}"

# Run the check
result = check_correctness()
print(result)
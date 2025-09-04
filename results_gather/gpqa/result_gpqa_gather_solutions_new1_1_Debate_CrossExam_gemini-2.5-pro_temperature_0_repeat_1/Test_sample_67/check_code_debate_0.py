import sympy

def check_oscillation_frequency():
    """
    Symbolically calculates the oscillation frequency of the z-component of the magnetic moment
    for a spin-1/2 particle undergoing Larmor precession and checks it against the provided answer.
    """
    # 1. Define symbols for the physical quantities
    gamma, B, t, hbar = sympy.symbols('gamma B t hbar', real=True, positive=True)

    # 2. Define the spin-1/2 operators using Pauli matrices
    # S_k = (hbar/2) * sigma_k
    sigma_y = sympy.Matrix([[0, -sympy.I], [sympy.I, 0]])
    sigma_z = sympy.Matrix([[1, 0], [0, -1]])
    S_y = (hbar / 2) * sigma_y
    S_z = (hbar / 2) * sigma_z

    # 3. Define the initial state: spin-up along the Z-axis
    # |psi(0)> = |+z> = [1, 0]^T
    psi_0 = sympy.Matrix([[1], [0]])

    # 4. Define the new Hamiltonian for a B-field in the +Y direction
    # H = -mu . B = -gamma * S . B = -gamma * B * S_y
    H = -gamma * B * S_y

    # 5. Calculate the time evolution operator U(t) = exp(-iHt/hbar)
    # This can be calculated directly using the matrix exponential.
    # U(t) = exp(-i * (-gamma * B * S_y) * t / hbar) = exp(i * gamma * B * t * S_y / hbar)
    # Let's use sympy's built-in matrix exponential function for robustness.
    U_t = sympy.exp(-sympy.I * H * t / hbar)
    U_t = sympy.simplify(U_t)

    # 6. Calculate the time-evolved state |psi(t)> = U(t)|psi(0)>
    psi_t = U_t * psi_0

    # 7. Calculate the expectation value of S_z at time t
    # <S_z(t)> = <psi(t)| S_z |psi(t)> = psi_t^dagger * S_z * psi_t
    # .H in sympy gives the conjugate transpose (dagger)
    exp_Sz_t = (psi_t.H * S_z * psi_t)[0] # The result is a 1x1 matrix, so we extract the element
    exp_Sz_t = sympy.simplify(exp_Sz_t)

    # The expectation value of the magnetic moment is <mu_z(t)> = gamma * <S_z(t)>
    # This has the same oscillation frequency as <S_z(t)>.
    # The result is of the form A*cos(omega*t). We need to extract omega.

    # 8. Extract the angular frequency (omega) from the expression
    calculated_frequency = 0
    # The expression for exp_Sz_t will contain a cosine term.
    # We find the argument of that cosine term.
    cos_terms = exp_Sz_t.find(sympy.cos)
    if not cos_terms:
        return f"Checker Error: Could not find a cosine term in the expectation value {exp_Sz_t}."
    
    # The argument of the cosine function is omega*t
    argument = list(cos_terms)[0].args[0]
    
    # The frequency omega is the coefficient of t in the argument
    calculated_frequency = argument.coeff(t)

    # 9. Compare the calculated frequency with the given options and the final answer
    options = {
        "A": gamma * B / 2,
        "B": gamma * B / sympy.sqrt(2),
        "C": sympy.sqrt(2) * gamma * B,
        "D": gamma * B
    }
    
    llm_final_answer_key = "D"
    
    # Check if the calculated frequency matches the correct physical result
    correct_frequency = gamma * B
    if sympy.simplify(calculated_frequency - correct_frequency) != 0:
        return f"Incorrect. The checker's calculation yielded a frequency of {calculated_frequency}, which does not match the known correct Larmor frequency of {correct_frequency}."

    # Check if the LLM's chosen answer corresponds to the correct frequency
    if sympy.simplify(options[llm_final_answer_key] - correct_frequency) == 0:
        return "Correct"
    else:
        return f"Incorrect. The correct oscillation frequency is {correct_frequency}. This corresponds to option D. The provided answer was '{llm_final_answer_key}', but the analysis leading to it was correct."

# Run the check
result = check_oscillation_frequency()
print(result)
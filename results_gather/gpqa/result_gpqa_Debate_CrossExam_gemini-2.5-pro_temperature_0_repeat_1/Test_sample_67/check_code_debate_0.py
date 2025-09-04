import sympy

def check_spin_oscillation_frequency():
    """
    Symbolically calculates the oscillation frequency of <Sz(t)> for a spin-1/2 particle
    precessing around the Y-axis and checks it against the provided answer.
    """
    # 1. Define symbolic variables for the physical quantities.
    # gamma: gyromagnetic ratio, B: magnetic field strength, t: time, hbar: reduced Planck constant
    gamma, B, t, hbar = sympy.symbols('gamma B t hbar', real=True, positive=True)
    I = sympy.I  # Imaginary unit

    # 2. Define spin operators in the Sz basis using Pauli matrices.
    # sigma_y = [[0, -i], [i, 0]], sigma_z = [[1, 0], [0, -1]]
    Sy = (hbar / 2) * sympy.Matrix([[0, -I], [I, 0]])
    Sz = (hbar / 2) * sympy.Matrix([[1, 0], [0, -1]])

    # 3. Define the Hamiltonian for t > 0.
    # The magnetic field is B along the +Y direction. H = -μ·B = -γ(S·B) = -γB*Sy.
    H = -gamma * B * Sy

    # 4. Define the initial state at t=0.
    # The particle is initially aligned with a field in the +Z direction,
    # so it is in the spin-up state |+Z>. In the Sz basis, this is [1, 0].
    psi_0 = sympy.Matrix([1, 0])

    # 5. Calculate the time evolution operator U(t) = exp(-iHt/ħ).
    # sympy.exp() computes the matrix exponential.
    U_t = sympy.exp(-I * H * t / hbar)
    U_t = sympy.simplify(U_t)

    # 6. Calculate the state of the system at time t: |ψ(t)> = U(t)|ψ(0)>.
    psi_t = U_t * psi_0

    # 7. Calculate the expectation value of Sz at time t: <Sz(t)> = <ψ(t)|Sz|ψ(t)>.
    # .H in sympy performs the conjugate transpose (Hermitian conjugate).
    expectation_Sz_t = psi_t.H * Sz * psi_t

    # The result is a 1x1 matrix; we extract the scalar value.
    expectation_Sz_t_scalar = sympy.simplify(expectation_Sz_t[0, 0])

    # The expected result from the derivation is (ħ/2) * cos(γBt). Let's verify.
    expected_Sz_expression = (hbar / 2) * sympy.cos(gamma * B * t)
    if sympy.simplify(expectation_Sz_t_scalar - expected_Sz_expression) != 0:
        return (f"The symbolic calculation of <Sz(t)> yielded an unexpected result. "
                f"Calculated: {expectation_Sz_t_scalar}, "
                f"Expected from derivation: {expected_Sz_expression}")

    # 8. Extract the angular frequency (ω) from the expression A*cos(ωt).
    # The argument of the cosine function is ω*t.
    # We can find ω by taking the derivative of the argument with respect to t.
    try:
        # Find the cosine term in the expression
        cos_term = [arg for arg in expectation_Sz_t_scalar.args if isinstance(arg, sympy.cos)][0]
        # Get its argument
        cos_argument = cos_term.args[0]
        # The frequency is the coefficient of t in the argument
        calculated_frequency = cos_argument.diff(t)
    except (IndexError, AttributeError):
        return f"Could not extract frequency from the expression: {expectation_Sz_t_scalar}"

    # The question asks for the oscillation frequency. The answer A is gamma*B.
    answer_frequency = gamma * B

    # 9. Compare the calculated frequency with the answer.
    if calculated_frequency == answer_frequency:
        # The average of the z-component of the moment is <μz> = γ<Sz>.
        # Since <Sz> oscillates with frequency γB, <μz> also oscillates with frequency γB.
        # The answer A is correct.
        return "Correct"
    else:
        return (f"The calculated oscillation frequency is {calculated_frequency}, "
                f"but the frequency from answer A is {answer_frequency}. The answer is incorrect.")

# Run the check
result = check_spin_oscillation_frequency()
print(result)
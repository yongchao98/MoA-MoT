import sympy

def check_larmor_frequency():
    """
    Symbolically calculates the oscillation frequency of the z-component of the spin
    for a spin-1/2 particle undergoing Larmor precession.

    The problem states:
    1. A spin-1/2 particle starts in a state aligned with a +Z magnetic field.
       This means its initial state is |+z>, the spin-up eigenstate of Sz.
    2. The field is suddenly switched to the +Y direction with the same magnitude B.
       The new Hamiltonian is H = -γB*Sy.
    3. We need to find the oscillation frequency of the expectation value <μz(t)>,
       which is proportional to <Sz(t)>.

    The code will perform this calculation and check if the result matches the
    provided answer.
    """
    try:
        # 1. Define symbolic variables for the physical quantities.
        gamma, B, t, hbar = sympy.symbols('gamma B t hbar', real=True, positive=True)
        I = sympy.I  # Imaginary unit

        # 2. Define operators and the initial state in the Sz basis.
        # Sy = (ħ/2) * PauliY
        Sy = (hbar / 2) * sympy.Matrix([[0, -I], [I, 0]])
        # Sz = (ħ/2) * PauliZ
        Sz = (hbar / 2) * sympy.Matrix([[1, 0], [0, -1]])
        # Initial state |ψ(0)> = |+z>
        psi_0 = sympy.Matrix([[1], [0]])

        # 3. Define the Hamiltonian for t > 0.
        H = -gamma * B * Sy

        # 4. Calculate the time evolution operator U = exp(-iHt/ħ).
        # sympy.exp can compute the exponential of a matrix.
        U = sympy.exp(-I * H * t / hbar)
        U = sympy.simplify(U)

        # 5. Calculate the time-evolved state |ψ(t)> = U|ψ(0)>.
        psi_t = U * psi_0

        # 6. Calculate the expectation value <Sz>(t) = <ψ(t)|Sz|ψ(t)>.
        # .H gives the Hermitian conjugate (conjugate transpose).
        exp_Sz_t = psi_t.H * Sz * psi_t
        
        # The result is a 1x1 matrix; extract the scalar value and simplify.
        exp_Sz_t_scalar = sympy.simplify(exp_Sz_t[0, 0])

        # 7. Extract the oscillation frequency from the result.
        # The result will be of the form A*cos(ω*t). We need to find ω.
        
        # Find the cosine function in the expression
        cos_term = None
        for arg in sympy.preorder_traversal(exp_Sz_t_scalar):
            if isinstance(arg, sympy.cos):
                cos_term = arg
                break
        
        if cos_term is None:
            return f"Calculation Error: The expectation value <Sz>(t) is not in the expected cosine form. Result: {exp_Sz_t_scalar}"

        # The argument of the cosine is ω*t
        cos_argument = cos_term.args[0]
        
        # The angular frequency ω is the derivative of the argument with respect to time t.
        calculated_frequency = sympy.diff(cos_argument, t)

        # 8. Compare the calculated frequency with the given answer.
        # The final answer provided is D, which corresponds to gamma*B.
        expected_frequency = gamma * B

        if calculated_frequency == expected_frequency:
            return "Correct"
        else:
            return f"Incorrect. The symbolic calculation shows the oscillation frequency is {calculated_frequency}, but the answer D corresponds to {expected_frequency}."

    except Exception as e:
        return f"An error occurred during the calculation: {e}"

# Run the check
result = check_larmor_frequency()
print(result)
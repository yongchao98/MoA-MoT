import sympy
from sympy import I, cos, sin, eye, simplify, diff

def check_correctness():
    """
    Checks the correctness of the LLM's answer for the spin precession problem.

    This function symbolically calculates the oscillation frequency of the expectation
    value of the z-component of a spin-1/2 particle's magnetic moment when the
    magnetic field is suddenly switched from the Z to the Y direction.
    """
    try:
        # 1. Define the physical symbols
        gamma, B, t, hbar = sympy.symbols('gamma B t hbar', real=True, positive=True)

        # 2. Define the spin operators and the initial state in the S_z basis.
        # The initial state is spin-up along Z, |+z> = [1, 0]^T.
        S_y = (hbar / 2) * sympy.Matrix([[0, -I], [I, 0]])
        S_z = (hbar / 2) * sympy.Matrix([[1, 0], [0, -1]])
        psi_0 = sympy.Matrix([1, 0])

        # 3. Define the new Hamiltonian for a B field in the +Y direction.
        H = -gamma * B * S_y

        # 4. Calculate the time evolution operator U(t) = exp(-iHt/ħ).
        # We use the identity exp(i*a*σ_k) = cos(a)*I + i*sin(a)*σ_k,
        # where σ_k is a Pauli matrix and I is the identity matrix.
        # The argument of our exponential is (i*gamma*B*t/ħ) * S_y = i*(gamma*B*t/2)*σ_y.
        # So, a = gamma*B*t/2 and k=y.
        theta = gamma * B * t / 2
        sigma_y = sympy.Matrix([[0, -I], [I, 0]])
        identity_matrix = eye(2)
        U = cos(theta) * identity_matrix + I * sin(theta) * sigma_y

        # 5. Calculate the time-evolved state |ψ(t)> = U(t)|ψ(0)>.
        psi_t = U * psi_0

        # 6. Calculate the expectation value <S_z(t)> = <ψ(t)|S_z|ψ(t)>.
        # <ψ(t)| is the Hermitian conjugate (conjugate transpose) of |ψ(t)>.
        exp_Sz_t = (psi_t.H * S_z * psi_t)[0, 0] # The result is a 1x1 matrix.
        exp_Sz_t_simplified = simplify(exp_Sz_t)

        # The result should be of the form A*cos(ωt). We need to extract ω.
        # Find all cosine terms in the simplified expression.
        cos_terms = exp_Sz_t_simplified.find(cos)
        if not cos_terms:
            return (f"Failed to verify. The calculated expectation value <S_z(t)> = "
                    f"{exp_Sz_t_simplified} does not contain a cosine term, so the "
                    f"frequency cannot be determined.")

        # The argument of the cosine function contains the time dependence.
        cos_argument = list(cos_terms)[0].args[0]

        # The angular frequency ω is the coefficient of t in the argument.
        # A robust way to extract this is to differentiate the argument w.r.t. time.
        calculated_frequency = diff(cos_argument, t)

        # 7. Compare the calculated frequency with the provided answer.
        # The final answer from the LLM is <<<C>>>, which corresponds to gamma*B.
        llm_answer_symbolic = gamma * B

        if simplify(calculated_frequency - llm_answer_symbolic) == 0:
            return "Correct"
        else:
            return (f"Incorrect. The calculated oscillation frequency is {calculated_frequency}. "
                    f"The provided answer corresponds to {llm_answer_symbolic}.")

    except Exception as e:
        return f"An error occurred during the verification process: {e}"

# Run the check
result = check_correctness()
print(result)
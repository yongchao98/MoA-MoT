import sympy
from sympy.physics.quantum import Dagger

def check_oscillation_frequency():
    """
    Verifies the oscillation frequency of a spin-1/2 particle's magnetic moment
    by symbolically calculating its time evolution in the Schrödinger picture.
    """
    try:
        # 1. Define symbolic variables for the physical quantities.
        gamma, B, t, hbar = sympy.symbols('gamma B t hbar', real=True, positive=True)
        i = sympy.I

        # 2. Define the spin operators for a spin-1/2 particle in the standard Z-basis.
        # These are the Pauli matrices multiplied by hbar/2.
        sigma_y = sympy.Matrix([[0, -i], [i, 0]])
        S_y = (hbar / 2) * sigma_y
        S_z = (hbar / 2) * sympy.Matrix([[1, 0], [0, -1]])

        # 3. Define the initial state of the particle.
        # It is aligned with the +Z direction, so it's in the spin-up state |+z>.
        # In the Z-basis, this is represented by the vector [1, 0].
        psi_0 = sympy.Matrix([1, 0])

        # 4. Define the Hamiltonian for the new magnetic field in the +Y direction.
        # H = -μ · B = -γS · B = -γB*S_y
        H = -gamma * B * S_y

        # 5. Calculate the time evolution operator U(t) = exp(-iHt/ħ).
        # We use the identity: exp(iθσ_k) = cos(θ)I + i*sin(θ)σ_k
        # The argument of the exponential is: -iHt/ħ = i(γBt/ħ)S_y = i(γBt/2)σ_y
        # So, θ = γBt/2
        theta = gamma * B * t / 2
        identity_matrix = sympy.eye(2)
        U_t = sympy.cos(theta) * identity_matrix + i * sympy.sin(theta) * sigma_y

        # 6. Calculate the state of the particle at time t: |ψ(t)> = U(t)|ψ(0)>
        psi_t = U_t * psi_0

        # 7. Calculate the expectation value of the z-component of the spin: <S_z>(t) = <ψ(t)|S_z|ψ(t)>
        # <ψ(t)| is the Hermitian conjugate (conjugate transpose) of |ψ(t)>.
        bra_psi_t = Dagger(psi_t)
        
        # The result of the inner product is a 1x1 matrix, so we extract the single element.
        exp_Sz_t = (bra_psi_t * S_z * psi_t)[0]

        # 8. Simplify the resulting expression.
        exp_Sz_t_simplified = sympy.simplify(exp_Sz_t)
        
        # The result is expected to be of the form A*cos(ωt). We need to find ω.
        # The simplified expression should be (ħ/2) * cos(γBt).
        
        # 9. Extract the angular frequency (ω) from the expression.
        # The frequency is the coefficient of the time variable 't' inside the cosine function.
        
        # Find the cosine term in the expression
        cos_term = None
        if exp_Sz_t_simplified.is_Mul:
            for term in exp_Sz_t_simplified.args:
                if isinstance(term, sympy.cos):
                    cos_term = term
                    break
        
        if cos_term is None:
            return f"Error: Could not find a cosine term in the calculated expectation value: {exp_Sz_t_simplified}"

        # The argument of the cosine is the first (and only) element of its .args tuple
        cos_argument = cos_term.args[0]
        
        # The angular frequency is the derivative of the argument with respect to time.
        calculated_frequency = sympy.diff(cos_argument, t)

        # 10. Compare the calculated frequency with the frequency from the provided answer.
        # The question's options are:
        # A) sqrt(2)*gamma*B
        # B) gamma*B/2
        # C) gamma*B
        # D) gamma*B/sqrt(2)
        # The provided answer is <<<C>>>, which corresponds to the frequency gamma*B.
        answer_frequency = gamma * B
        
        # Check if the calculated frequency matches the answer's frequency.
        if sympy.simplify(calculated_frequency - answer_frequency) == 0:
            return "Correct"
        else:
            return (f"Incorrect: The calculated oscillation frequency is {calculated_frequency}, "
                    f"but the answer corresponds to a frequency of {answer_frequency}.")

    except ImportError:
        return "Skipping check: The 'sympy' library is required but not installed."
    except Exception as e:
        return f"An error occurred during the check: {e}"

# Run the check and print the result.
result = check_oscillation_frequency()
print(result)
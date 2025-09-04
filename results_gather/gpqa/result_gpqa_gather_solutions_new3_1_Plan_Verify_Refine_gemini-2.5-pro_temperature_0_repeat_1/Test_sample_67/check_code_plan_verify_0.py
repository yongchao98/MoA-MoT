import sympy

def check_larmor_frequency():
    """
    Symbolically calculates the oscillation frequency of <S_z(t)> for a spin-1/2 particle
    precessing in a magnetic field along the Y-axis.
    """
    try:
        # 1. Define symbolic variables
        gamma, B, t, hbar = sympy.symbols('gamma B t hbar', real=True, positive=True)
        I = sympy.I

        # 2. Set up operators and initial state
        # Pauli matrices
        sigma_y = sympy.Matrix([[0, -I], [I, 0]])
        sigma_z = sympy.Matrix([[1, 0], [0, -1]])

        # Spin operators
        S_y = (hbar / 2) * sigma_y
        S_z = (hbar / 2) * sigma_z

        # Initial state |psi(0)> = |+z> (spin-up along z-axis)
        psi_0 = sympy.Matrix([1, 0])

        # 3. Define the Hamiltonian for t > 0
        # H = -μ . B = -γS . (B * y_hat) = -γB * S_y
        H = -gamma * B * S_y

        # 4. Calculate the time evolution operator U(t) = exp(-iHt/ħ)
        # sympy.exp can compute the matrix exponential
        U = sympy.exp(H * t * (-I / hbar))

        # 5. Calculate the state at time t: |ψ(t)> = U(t)|ψ(0)>
        psi_t = U * psi_0

        # 6. Calculate the expectation value of S_z at time t: <S_z(t)> = <ψ(t)|S_z|ψ(t)>
        # .H gives the conjugate transpose (Hermitian conjugate)
        exp_Sz_matrix = psi_t.H * S_z * psi_t
        
        # The result is a 1x1 matrix, extract the scalar value
        exp_Sz = sympy.simplify(exp_Sz_matrix[0, 0])

        # 7. Extract the oscillation frequency
        # The result is of the form A * cos(ωt). We find ω.
        
        # Find the cosine term in the expression
        cos_term = [arg for arg in exp_Sz.args if isinstance(arg, sympy.cos)]
        if not cos_term:
            return f"Calculation Error: The expectation value <S_z(t)> = {exp_Sz} is not in the expected oscillatory form."
        
        # Get the argument of the cosine function, e.g., gamma*B*t
        argument = cos_term[0].args[0]
        
        # The angular frequency ω is the coefficient of t, which can be found by differentiating the argument w.r.t. t
        calculated_frequency = sympy.diff(argument, t)

        # 8. Compare with the provided answer
        # The provided answer is D, which corresponds to the expression gamma*B
        expected_frequency = gamma * B
        
        if sympy.simplify(calculated_frequency - expected_frequency) == 0:
            return "Correct"
        else:
            return (f"Incorrect. The calculated oscillation frequency is {calculated_frequency}, "
                    f"but the frequency corresponding to answer D is {expected_frequency}.")

    except Exception as e:
        return f"An error occurred during the calculation: {e}"

# Run the check
result = check_larmor_frequency()
print(result)
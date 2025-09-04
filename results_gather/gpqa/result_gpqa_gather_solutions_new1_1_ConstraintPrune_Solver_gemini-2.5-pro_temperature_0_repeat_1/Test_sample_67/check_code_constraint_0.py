import sympy
from sympy import I, cos, sin, Matrix, symbols, simplify, diff

def check_correctness():
    """
    This function symbolically derives the oscillation frequency for the given quantum mechanics problem
    and checks if it matches the frequency corresponding to the provided answer 'A'.

    The problem involves a spin-1/2 particle, initially in a spin-up state along the Z-axis,
    suddenly subjected to a magnetic field in the Y-direction. We need to find the oscillation
    frequency of the expectation value of the z-component of its magnetic moment, <μ_z>.
    """
    try:
        # 1. Define symbols for physical quantities
        gamma, B, t, hbar = symbols('gamma B t hbar', real=True, positive=True)

        # 2. Define spin operators and the initial state in the S_z basis
        # Pauli matrices are the core of spin-1/2 operators
        sigma_y = Matrix([[0, -I], [I, 0]])
        sigma_z = Matrix([[1, 0], [0, -1]])
        
        # Spin operators are proportional to Pauli matrices
        Sy = (hbar / 2) * sigma_y
        Sz = (hbar / 2) * sigma_z
        
        # The initial state is spin-up along Z, represented by the vector [1, 0]
        psi_0 = Matrix([1, 0])

        # 3. Construct the Hamiltonian for the new magnetic field
        # H = -μ · B = -γS · B = -γB S_y
        H = -gamma * B * Sy

        # 4. Calculate the time-evolved state |ψ(t)> = U(t)|ψ(0)>
        # The time evolution operator is U(t) = exp(-iHt/ħ).
        # Using the identity exp(i*θ*σ_k) = cos(θ)*I + i*sin(θ)*σ_k, where θ = γBt/2
        theta = gamma * B * t / 2
        identity_matrix = Matrix([[1, 0], [0, 1]])
        U_t = cos(theta) * identity_matrix + I * sin(theta) * sigma_y
        psi_t = U_t * psi_0

        # 5. Calculate the expectation value <S_z(t)> = <ψ(t)|S_z|ψ(t)>
        # The .H attribute gives the Hermitian conjugate (conjugate transpose)
        expectation_Sz_t = (psi_t.H * Sz * psi_t)[0]  # The result is a 1x1 matrix
        expectation_Sz_t = simplify(expectation_Sz_t)

        # The expected result is of the form A*cos(ωt). We need to extract ω.
        # The argument of the cosine function is ωt. The frequency ω is the derivative of the argument with respect to t.
        
        # Find the cosine term in the simplified expression
        cos_term = None
        if isinstance(expectation_Sz_t, sympy.cos):
            cos_term = expectation_Sz_t
        elif isinstance(expectation_Sz_t, sympy.mul.Mul):
            for arg in expectation_Sz_t.args:
                if isinstance(arg, sympy.cos):
                    cos_term = arg
                    break
        
        if cos_term is None:
            return f"Could not find the oscillating term in the expectation value. Derived result was: {expectation_Sz_t}"

        argument = cos_term.args[0]
        calculated_frequency = diff(argument, t)

        # 6. Compare the calculated frequency with the frequency from option A
        # Option A is gamma*B
        frequency_A = gamma * B

        if simplify(calculated_frequency - frequency_A) == 0:
            # The derived frequency matches the frequency from option A.
            # The oscillation frequency of <μ_z> is the same as for <S_z>.
            return "Correct"
        else:
            return f"The derived oscillation frequency is {calculated_frequency}, which does not match the frequency from option A ({frequency_A})."

    except Exception as e:
        return f"An error occurred during the symbolic calculation: {e}"

# Run the check
result = check_correctness()
print(result)
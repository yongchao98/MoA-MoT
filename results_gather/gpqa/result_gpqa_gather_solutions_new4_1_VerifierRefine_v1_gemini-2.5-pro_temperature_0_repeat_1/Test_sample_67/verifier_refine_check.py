import sympy
from sympy import I, Matrix, cos, sin, symbols

def check_correctness():
    """
    Symbolically calculates the oscillation frequency of the z-component of the magnetic moment
    for a spin-1/2 particle undergoing Larmor precession.
    """
    try:
        # 1. Define symbols for physical quantities
        gamma, B, t, hbar = symbols('gamma B t hbar', real=True, positive=True)

        # 2. Define Pauli matrices and spin operators
        # We use the standard basis where |+z> = [1, 0]^T and |-z> = [0, 1]^T
        sigma_y = Matrix([[0, -I], [I, 0]])
        sigma_z = Matrix([[1, 0], [0, -1]])
        
        Sy = (hbar / 2) * sigma_y
        Sz = (hbar / 2) * sigma_z

        # 3. Define the initial state |+z>
        psi_0 = Matrix([1, 0])

        # 4. Define the new Hamiltonian H = -gamma * B * Sy
        H = -gamma * B * Sy

        # 5. Calculate the time evolution of the state |ψ(t)>
        # The state evolves as |ψ(t)> = exp(-iHt/ħ)|ψ(0)>
        # The operator is exp(-i*(-γBSy)t/ħ) = exp(iγBtSy/ħ)
        # Using the identity exp(iθσ_k) = cos(θ)I + i*sin(θ)σ_k, with θ = γBt/2
        theta = gamma * B * t / 2
        identity_matrix = Matrix([[1, 0], [0, 1]])
        U = cos(theta) * identity_matrix + I * sin(theta) * sigma_y
        psi_t = U * psi_0

        # 6. Calculate the expectation value of Sz at time t
        # <Sz>(t) = <ψ(t)|Sz|ψ(t)> = ψ(t)† * Sz * ψ(t)
        exp_Sz_t = (psi_t.H * Sz * psi_t)[0, 0] # Result is a 1x1 matrix
        exp_Sz_t_simplified = sympy.simplify(exp_Sz_t)

        # 7. The expectation value of the magnetic moment is γ * <Sz>(t)
        exp_mu_z_t = gamma * exp_Sz_t_simplified

        # Expected form is A*cos(ωt). Let's verify.
        # The simplified result should be (γħ/2) * cos(γBt)
        expected_result = (gamma * hbar / 2) * cos(gamma * B * t)
        if sympy.simplify(exp_mu_z_t - expected_result) != 0:
            return f"Calculation error: The derived expectation value <μz>(t) = {exp_mu_z_t} does not match the expected form."

        # 8. Extract the oscillation frequency (ω) from the argument of the cosine
        # The argument is γBt. The frequency is the coefficient of t.
        cos_term = [arg for arg in exp_mu_z_t.args if isinstance(arg, sympy.cos)][0]
        cos_argument = cos_term.args[0]
        calculated_frequency = sympy.diff(cos_argument, t)

        # 9. Compare with the frequency from the proposed answer (C: γB)
        proposed_frequency = gamma * B

        if calculated_frequency == proposed_frequency:
            return "Correct"
        else:
            return f"The calculated oscillation frequency is {calculated_frequency}, but the answer 'C' implies a frequency of {proposed_frequency}."

    except Exception as e:
        return f"An error occurred during the symbolic calculation: {e}"

# Run the check
result = check_correctness()
print(result)
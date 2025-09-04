import sympy

def check_larmor_precession_frequency():
    """
    Symbolically calculates the oscillation frequency of <Sz> for a spin-1/2 particle
    precessing around the Y-axis.
    """
    try:
        # 1. Define symbols for the physical quantities
        gamma, B, t, hbar = sympy.symbols('gamma B t hbar', real=True, positive=True)
        I = sympy.I

        # 2. Define the spin-1/2 operators using Pauli matrices
        # Sy = (hbar/2) * sigma_y
        # Sz = (hbar/2) * sigma_z
        S_y = (hbar / 2) * sympy.Matrix([[0, -I], [I, 0]])
        S_z = (hbar / 2) * sympy.Matrix([[1, 0], [0, -1]])

        # 3. Define the initial state: spin-up along Z, |+z>
        psi_initial = sympy.Matrix([[1], [0]])

        # 4. Define the new Hamiltonian for a B-field in the +Y direction
        H = -gamma * B * S_y

        # 5. Calculate the time evolution operator U(t) = exp(-iHt/hbar)
        # The argument of the matrix exponential is M = -I * H * t / hbar
        M = -I * H * t / hbar
        # U(t) is the matrix exponential of M
        U = sympy.exp(M)
        U.simplify()

        # 6. Calculate the state at time t: |psi(t)> = U(t)|psi(0)>
        psi_t = U * psi_initial
        psi_t.simplify()

        # 7. Calculate the expectation value <Sz(t)> = <psi(t)|Sz|psi(t)>
        # <psi(t)| is the conjugate transpose of |psi(t)>
        bra_psi_t = psi_t.H
        expectation_Sz_matrix = bra_psi_t * S_z * psi_t
        
        # The result is a 1x1 matrix, so we extract the scalar element
        expectation_Sz = expectation_Sz_matrix[0, 0].simplify()

        # 8. Extract the oscillation frequency from the result
        # The result should be of the form A*cos(omega*t).
        # We find the argument of the cosine function.
        
        # Find the cosine term in the expression
        cos_term = [arg for arg in sympy.preorder_traversal(expectation_Sz) if isinstance(arg, sympy.cos)]
        
        if not cos_term:
            return f"Calculation failed to produce an oscillating term. Result: {expectation_Sz}"
            
        # Get the argument of the cosine function
        cos_argument = cos_term[0].args[0]
        
        # The angular frequency is the derivative of the argument with respect to time
        calculated_frequency = sympy.diff(cos_argument, t)

        # The answer given is D, which corresponds to the frequency gamma*B
        answer_frequency = gamma * B

        # 9. Check if the calculated frequency matches the answer
        if calculated_frequency == answer_frequency:
            # The oscillation frequency of the magnetic moment <mu_z> is the same as for <S_z>
            # because mu_z = gamma * S_z, and gamma is a constant.
            # The provided answer D corresponds to gamma*B, which matches our calculation.
            return "Correct"
        else:
            return (f"The calculated oscillation frequency is {calculated_frequency}, "
                    f"but the answer 'D' implies a frequency of {answer_frequency}.")

    except Exception as e:
        return f"An error occurred during the calculation: {e}"

# Run the check
result = check_larmor_precession_frequency()
print(result)
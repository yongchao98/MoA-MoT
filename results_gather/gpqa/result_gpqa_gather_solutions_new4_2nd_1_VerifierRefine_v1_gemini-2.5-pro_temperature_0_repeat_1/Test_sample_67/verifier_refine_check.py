import sympy

def check_correctness():
    """
    This function symbolically calculates the oscillation frequency of the z-component
    of a spin-1/2 particle's magnetic moment in a magnetic field along the Y-axis.
    It verifies if the frequency is gamma*B, as stated in the provided answer.
    """
    try:
        # 1. Define symbols for the physical quantities
        gamma, B, t, hbar = sympy.symbols('gamma B t hbar', real=True, positive=True)
        I = sympy.I  # Imaginary unit

        # 2. Define spin-1/2 operators using Pauli matrices
        # Sz = (hbar/2) * sigma_z
        # Sy = (hbar/2) * sigma_y
        sigma_y = sympy.Matrix([[0, -I], [I, 0]])
        sigma_z = sympy.Matrix([[1, 0], [0, -1]])
        
        Sy = hbar / 2 * sigma_y
        Sz = hbar / 2 * sigma_z

        # 3. Define the Hamiltonian for t > 0 (field along +Y)
        H = -gamma * B * Sy

        # 4. Define the initial state at t=0
        # The particle is aligned with the Z-field, so it's in the |+z> state.
        psi_0 = sympy.Matrix([1, 0])

        # 5. Calculate the time evolution operator U(t) = exp(-iHt/hbar)
        # For a spin-1/2 system, U(t) can be calculated using the identity:
        # exp(i*theta*sigma_k) = cos(theta)*I + i*sin(theta)*sigma_k
        # The argument of the exponential is: -I * H * t / hbar
        # = -I * (-gamma * B * Sy) * t / hbar
        # = I * gamma * B * t * (hbar/2 * sigma_y) / hbar
        # = I * (gamma * B * t / 2) * sigma_y
        theta = gamma * B * t / 2
        identity_matrix = sympy.eye(2)
        U = sympy.cos(theta) * identity_matrix + I * sympy.sin(theta) * sigma_y

        # 6. Calculate the state at time t: |psi(t)> = U(t)|psi(0)>
        psi_t = U * psi_0

        # 7. Calculate the expectation value of Sz at time t: <Sz>(t) = <psi(t)|Sz|psi(t)>
        # <psi(t)| is the conjugate transpose of |psi(t)>
        exp_Sz_t = (psi_t.H * Sz * psi_t)[0]  # .H is conjugate transpose, result is a 1x1 matrix

        # 8. Simplify the expression
        simplified_exp_Sz_t = sympy.simplify(exp_Sz_t)
        
        # The expected result is of the form A*cos(omega*t)
        # Let's find the frequency 'omega' from the expression
        
        # Find the cosine term in the simplified expression
        cos_term = None
        if isinstance(simplified_exp_Sz_t, sympy.cos):
             cos_term = simplified_exp_Sz_t
        elif isinstance(simplified_exp_Sz_t, sympy.Mul):
            for arg in simplified_exp_Sz_t.args:
                if isinstance(arg, sympy.cos):
                    cos_term = arg
                    break
        
        if cos_term is None:
            return f"Calculation Error: The expectation value <Sz>(t) = {simplified_exp_Sz_t} is not in the expected oscillatory form."

        # The argument of the cosine function is the phase, phi(t)
        phase = cos_term.args[0]
        
        # The angular frequency is the derivative of the phase with respect to time
        calculated_frequency = sympy.diff(phase, t)
        
        # 9. Compare with the expected frequency from the answer
        # The answer states the frequency is gamma*B (Option C)
        expected_frequency = gamma * B

        if sympy.simplify(calculated_frequency - expected_frequency) == 0:
            return "Correct"
        else:
            return f"Incorrect: The calculated oscillation frequency is {calculated_frequency}, but the answer claims it is {expected_frequency}."

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Run the check
result = check_correctness()
print(result)
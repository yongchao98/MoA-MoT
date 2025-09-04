import sympy

def check_spin_oscillation_frequency():
    """
    Symbolically calculates the oscillation frequency of the z-component of the magnetic moment
    to verify the correctness of the given answer.
    """
    # 1. Define symbolic variables for the physical quantities.
    # gamma: gyromagnetic ratio, B: magnetic field magnitude, hbar: reduced Planck constant, t: time
    gamma, B, hbar, t = sympy.symbols('gamma B hbar t', real=True, positive=True)

    # 2. Define the spin-1/2 operators using Pauli matrices.
    # S_k = (hbar/2) * sigma_k
    sigma_y = sympy.Matrix([[0, -sympy.I], [sympy.I, 0]])
    sigma_z = sympy.Matrix([[1, 0], [0, -1]])
    
    Sy = (hbar / 2) * sigma_y
    Sz = (hbar / 2) * sigma_z

    # 3. Define the initial state at t=0.
    # The particle is aligned with a B field in the +Z direction, so it's in the |+z> eigenstate.
    psi_0 = sympy.Matrix([1, 0])

    # 4. Define the Hamiltonian for t > 0.
    # The magnetic field is B in the +Y direction. H = -μ ⋅ B = - (gamma * S) ⋅ (B * y_hat)
    H = -gamma * B * Sy

    # 5. Calculate the time evolution operator U(t) = exp(-i*H*t/hbar).
    # The argument of the exponential is (-i * (-gamma * B * Sy) * t / hbar) = i * gamma * B * t * Sy / hbar
    # Let omega = gamma * B. The argument becomes i * omega * t * Sy / hbar.
    # Since Sy = (hbar/2)*sigma_y, the argument simplifies to i * (omega * t / 2) * sigma_y.
    # We use the identity exp(i*a*sigma_k) = cos(a)*I + i*sin(a)*sigma_k
    omega = gamma * B
    angle = omega * t / 2
    U_t = sympy.cos(angle) * sympy.eye(2) + sympy.I * sympy.sin(angle) * sigma_y

    # 6. Calculate the state vector at time t: |ψ(t)> = U(t)|ψ(0)>
    psi_t = U_t * psi_0

    # 7. Calculate the expectation value of Sz at time t: <Sz>(t) = <ψ(t)|Sz|ψ(t)>
    # In sympy, .H performs the conjugate transpose (Hermitian conjugate).
    sz_expectation_t = sympy.simplify((psi_t.H * Sz * psi_t)[0, 0])

    # The result should be (hbar/2) * cos(omega*t). Let's verify.
    expected_form = (hbar / 2) * sympy.cos(omega * t)
    if not sympy.simplify(sz_expectation_t - expected_form) == 0:
        return f"Calculation error: The derived expectation value {sz_expectation_t} does not match the expected form {expected_form}."

    # 8. Extract the angular frequency from the result.
    # The expression for <Sz>(t) is A*cos(freq*t). The angular frequency is 'freq'.
    # We can get this by finding the argument of the cosine function and dividing by t.
    if not isinstance(sz_expectation_t, sympy.Mul) or not isinstance(sz_expectation_t.args[-1], sympy.cos):
         return f"The expectation value {sz_expectation_t} is not in the expected A*cos(omega*t) form."
    
    cos_arg = sz_expectation_t.args[-1].args[0]
    calculated_frequency = cos_arg / t

    # 9. Compare the calculated frequency with the provided answer (Option C).
    # Option C is gamma*B, which we have defined as omega.
    answer_C_frequency = gamma * B

    if sympy.simplify(calculated_frequency - answer_C_frequency) == 0:
        return "Correct"
    else:
        return (f"The derived oscillation frequency is {calculated_frequency}, but the answer C "
                f"corresponds to a frequency of {answer_C_frequency}. The answer is incorrect.")

# Run the check
result = check_spin_oscillation_frequency()
print(result)
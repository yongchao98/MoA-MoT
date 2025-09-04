import sympy
from sympy import I, sin, cos, exp

def check_correctness():
    """
    Verifies the oscillation frequency of a spin-1/2 particle's magnetic moment.

    The code solves the time-dependent Schrödinger equation for the given scenario
    to find the expectation value of the z-component of the magnetic moment, <μ_z>(t).
    It then extracts the angular frequency from this time-dependent expression and
    compares it to the frequency given in the answer.

    Returns:
        str: "Correct" if the answer is correct, otherwise a string explaining the discrepancy.
    """
    # 1. Define the physical constants and variables as symbolic objects.
    # We can set hbar=1 as it cancels out in the final frequency calculation,
    # but we'll keep it for completeness.
    gamma, B, t, hbar = sympy.symbols('gamma B t hbar', real=True, positive=True)

    # 2. Define the spin operators in the standard basis (eigenbasis of Sz).
    # Sz = (hbar/2) * sigma_z, Sy = (hbar/2) * sigma_y
    # We only need Sy for the Hamiltonian and Sz for the expectation value.
    sigma_y = sympy.Matrix([[0, -I], [I, 0]])
    sigma_z = sympy.Matrix([[1, 0], [0, -1]])
    
    Sy = (hbar / 2) * sigma_y
    Sz = (hbar / 2) * sigma_z
    mu_z = gamma * Sz

    # 3. Define the initial state.
    # The particle is initially aligned with the +Z direction, so its state is
    # the "spin up" eigenstate of Sz, which is |+z> = [1, 0]^T.
    psi_0 = sympy.Matrix([1, 0])

    # 4. Define the Hamiltonian for the new magnetic field along the +Y direction.
    # H = -μ . B = - (gamma * S) . (B * y_hat) = -gamma * B * Sy
    H = -gamma * B * Sy

    # 5. Calculate the time evolution operator U(t) = exp(-iHt/hbar).
    # U(t) = exp(-I * (-gamma * B * Sy) * t / hbar)
    # U(t) = exp(I * gamma * B * t * Sy / hbar)
    # Substituting Sy = (hbar/2)*sigma_y, the hbar cancels:
    # U(t) = exp(I * (gamma * B * t / 2) * sigma_y)
    # Sympy can compute the matrix exponential directly.
    U_t = sympy.exp(I * gamma * B * t * Sy / hbar)

    # 6. Calculate the state of the system at time t.
    # |ψ(t)> = U(t)|ψ(0)>
    psi_t = U_t * psi_0

    # 7. Calculate the expectation value of μ_z at time t.
    # <μ_z>(t) = <ψ(t)|μ_z|ψ(t)>
    # <ψ(t)| is the conjugate transpose of |ψ(t)>, denoted by .H in sympy.
    # The result is a 1x1 matrix, so we extract the single element.
    exp_mu_z_t = (psi_t.H * mu_z * psi_t)[0]

    # 8. Simplify the expression to identify the oscillation.
    # The expression simplifies to a cosine function.
    # cos(A)^2 - sin(A)^2 = cos(2A)
    exp_mu_z_t_simplified = sympy.simplify(exp_mu_z_t)
    
    # The simplified expression is (gamma*hbar/2)*cos(gamma*B*t).
    # This shows an oscillation of the form cos(ωt).

    # 9. Extract the angular frequency (ω) from the argument of the cosine.
    # The argument of the cosine function is ω*t.
    try:
        # Find all functions in the expression
        trig_functions = exp_mu_z_t_simplified.atoms(sympy.Function)
        # Filter for cosine or sine
        cos_or_sin_func = next(f for f in trig_functions if isinstance(f, (sympy.cos, sympy.sin)))
        # The argument of the function is the first element in its args tuple
        argument = cos_or_sin_func.args[0]
        # The angular frequency is the coefficient of t in the argument
        calculated_omega = argument.coeff(t)
    except (StopIteration, AttributeError):
        return "Error: Could not extract the oscillation frequency from the calculated expectation value."

    # 10. Compare the calculated frequency with the proposed answer.
    # The proposed answer is B) gamma*B
    answer_omega = gamma * B

    if sympy.simplify(calculated_omega - answer_omega) == 0:
        return "Correct"
    else:
        return (f"The answer is incorrect. The calculated oscillation frequency is {calculated_omega}, "
                f"but the proposed answer corresponds to a frequency of {answer_omega}.")

# Run the check
result = check_correctness()
print(result)
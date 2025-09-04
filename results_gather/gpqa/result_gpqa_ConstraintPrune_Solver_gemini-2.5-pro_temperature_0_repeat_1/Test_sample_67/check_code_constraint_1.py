import sympy
from sympy import I, sin, cos, exp, trigsimp

def check_correctness():
    """
    Checks the correctness of the given answer for the quantum mechanics problem.

    The problem asks for the oscillation frequency of the z-component of a spin-1/2
    particle's magnetic moment after a magnetic field is switched from the +Z to
    the +Y direction.
    """
    try:
        # 1. Define symbolic variables for the physical quantities.
        # gamma: gyromagnetic ratio, B: magnetic field strength, t: time, hbar: reduced Planck constant.
        gamma, B, t, hbar = sympy.symbols('gamma B t hbar', real=True, positive=True)

        # The system is a spin-1/2 particle, described by 2x2 matrices in the |+z>, |-z> basis.
        # Define the necessary spin operators using Pauli matrices (S = hbar/2 * sigma).
        sigma_y = sympy.Matrix([[0, -I], [I, 0]])
        sigma_z = sympy.Matrix([[1, 0], [0, -1]])
        S_y = (hbar / 2) * sigma_y
        S_z = (hbar / 2) * sigma_z

        # 2. Set the initial state.
        # The particle is initially aligned with a field in the +Z direction,
        # so its state is the spin-up eigenstate of Sz, |+z>.
        psi_0 = sympy.Matrix([1, 0])

        # 3. Define the Hamiltonian for t > 0.
        # The magnetic field is along +Y, so H = -mu . B = -gamma * S . B = -gamma * B * S_y.
        H = -gamma * B * S_y

        # 4. Calculate the time evolution of the state.
        # The time evolution operator is U(t) = exp(-iHt/hbar).
        # U(t) = exp(-I * (-gamma * B * S_y) * t / hbar) = exp(I * gamma * B * t * S_y / hbar)
        # sympy.exp() can compute the matrix exponential.
        U_t = sympy.exp(I * gamma * B * t * S_y / hbar)
        
        # The time-evolved state is |psi(t)> = U(t)|psi(0)>.
        psi_t = trigsimp(U_t * psi_0)

        # 5. Calculate the expectation value of the z-component of the magnetic moment.
        # <mu_z(t)> = <psi(t)| mu_z |psi(t)> = gamma * <psi(t)| Sz |psi(t)>
        # <psi(t)| is the conjugate transpose of |psi(t)>, which is psi_t.H in sympy.
        exp_val_Sz_matrix = psi_t.H * S_z * psi_t
        
        # The result is a 1x1 matrix; extract the scalar value and simplify.
        exp_val_Sz = trigsimp(exp_val_Sz_matrix[0, 0])
        exp_val_mu_z = gamma * exp_val_Sz

        # The resulting expression is gamma * (hbar/2) * cos(gamma*B*t).
        # This shows an oscillation.

        # 6. Extract the oscillation frequency.
        # The expression is of the form A*cos(omega*t). The angular frequency is omega.
        # We can find omega by isolating the argument of the cosine function and
        # differentiating it with respect to time t.
        
        # Find all cosine terms in the expression.
        cos_terms = exp_val_mu_z.atoms(sympy.cos)
        if len(cos_terms) != 1:
            return f"The expectation value {exp_val_mu_z} is not a simple cosine function, so a single frequency cannot be determined."
        
        cos_term = list(cos_terms)[0]
        cos_argument = cos_term.args[0]
        
        # The angular frequency omega is d(argument)/dt.
        calculated_omega = sympy.diff(cos_argument, t)
        
        # 7. Compare with the given answer.
        # The answer A corresponds to a frequency of gamma*B.
        answer_omega = gamma * B
        
        if calculated_omega == answer_omega:
            return "Correct"
        else:
            return f"The calculated oscillation frequency is {calculated_omega}, but the answer A implies it is {answer_omega}. Therefore, the answer is incorrect."

    except Exception as e:
        return f"An error occurred during the verification process: {e}"

# Run the check and print the result.
result = check_correctness()
if result == "Correct":
    print("Correct")
else:
    print(f"The answer is incorrect.\nReason: {result}")

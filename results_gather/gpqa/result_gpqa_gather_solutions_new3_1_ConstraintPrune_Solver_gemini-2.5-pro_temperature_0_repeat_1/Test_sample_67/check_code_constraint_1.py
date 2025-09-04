import sympy
from sympy import I, sin, cos, Matrix, Symbol

def check_larmor_frequency():
    """
    This function symbolically calculates the oscillation frequency of the z-component
    of a spin-1/2 particle's magnetic moment and verifies it against the provided answer.

    The problem setup is:
    1. Initial state: Spin-up along the Z-axis.
    2. New Hamiltonian: H = -gamma * B * S_y (due to B field in +Y direction).

    The function calculates <S_z(t)> and extracts the oscillation frequency.
    """
    try:
        # Define the physical symbols
        gamma = Symbol('gamma', real=True)
        B = Symbol('B', real=True, positive=True)
        t = Symbol('t', real=True, positive=True)
        hbar = Symbol('hbar', real=True, positive=True)

        # Define the Pauli Y matrix (sigma_y)
        sigma_y = Matrix([[0, -I], [I, 0]])
        
        # The spin-y operator is S_y = (hbar/2) * sigma_y
        # The Hamiltonian is H = -gamma * B * S_y
        # The time evolution operator is U(t) = exp(-i*H*t/hbar)
        # U(t) = exp(-i * (-gamma * B * S_y) * t / hbar)
        # U(t) = exp(i * gamma * B * t * S_y / hbar)
        # U(t) = exp(i * (gamma * B * t / 2) * sigma_y)
        
        # We use the identity: exp(i*theta*sigma_k) = cos(theta)*I + i*sin(theta)*sigma_k
        # where theta = gamma * B * t / 2
        theta = gamma * B * t / 2
        identity_matrix = Matrix([[1, 0], [0, 1]])
        U_t = cos(theta) * identity_matrix + I * sin(theta) * sigma_y

        # Initial state is spin-up along Z: |+z> = [1, 0]^T
        psi_0 = Matrix([1, 0])

        # The time-evolved state is |psi(t)> = U(t) * |psi(0)>
        psi_t = U_t * psi_0

        # The S_z operator is S_z = (hbar/2) * sigma_z
        sigma_z = Matrix([[1, 0], [0, -1]])
        S_z = (hbar / 2) * sigma_z

        # Calculate the expectation value of S_z at time t:
        # <S_z>(t) = <psi(t)| S_z |psi(t)> = psi(t).H * S_z * psi(t)
        # .H is the conjugate transpose (Hermitian conjugate)
        exp_Sz_t_matrix = psi_t.H * S_z * psi_t
        
        # The result is a 1x1 matrix, so we extract the single element and simplify
        exp_Sz_t = sympy.simplify(exp_Sz_t_matrix[0, 0])

        # The result should be of the form A*cos(omega*t).
        # From the calculation, exp_Sz_t = (hbar/2) * cos(gamma*B*t)
        # The argument of the cosine function is (gamma*B*t).
        # The angular frequency 'omega' is the coefficient of 't' in the argument.
        
        # We can extract the argument from the cosine function
        if not isinstance(exp_Sz_t.func, cos):
             return f"Calculation failed: The expectation value <S_z(t)> is not a simple cosine function. Result: {exp_Sz_t}"
        
        argument = exp_Sz_t.args[0]
        
        # The frequency is the argument divided by time 't'
        derived_frequency = sympy.simplify(argument / t)
        
        expected_frequency = gamma * B
        
        if sympy.simplify(derived_frequency - expected_frequency) != 0:
            return f"Derivation Error: The derived frequency is {derived_frequency}, but it should be {expected_frequency}."

        # The options given in the question
        options = {
            "A": gamma * B / sympy.sqrt(2),
            "B": gamma * B,
            "C": sympy.sqrt(2) * gamma * B,
            "D": gamma * B / 2
        }

        # The final answer provided by the LLM is 'B'
        provided_answer_key = 'B'
        
        # Check if our derived frequency matches the expression for the provided answer key
        if sympy.simplify(derived_frequency - options[provided_answer_key]) == 0:
            return "Correct"
        else:
            return (f"Incorrect. The derived oscillation frequency is {derived_frequency}. "
                    f"The given answer '{provided_answer_key}' corresponds to the expression {options[provided_answer_key]}, which is incorrect.")

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Run the check
result = check_larmor_frequency()
print(result)
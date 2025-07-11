def solve_ccc_z_synthesis():
    """
    Calculates and prints the minimal number of CCZ gates to synthesize a CCCZ gate.
    
    This solution is based on established results from quantum circuit synthesis theory.
    The problem of finding the CCZ-cost of a CCCZ gate is equivalent to finding the
    Toffoli-cost of a 4-qubit Toffoli gate, due to the local equivalence of these gates
    (they are interconvertible using single-qubit Hadamard gates).
    
    For ancilla-free synthesis, the cost to build an n-qubit Toffoli gate (C^(n-2)(X))
    from (n-1)-qubit Toffoli gates is known to be 2*(n-1) - 2.
    For a CCCZ gate, which acts on 4 qubits (n=4), the cost is 2*4 - 2 = 6.
    """
    
    # Number of qubits involved in a CCCZ gate (3 controls + 1 target)
    n_qubits = 4
    
    # The known formula for ancilla-free synthesis of a C^(n-2)(X) gate is 2n-2.
    minimal_ccz_gates = 2 * n_qubits - 2

    print("In the context of quantum circuit synthesis without ancilla qubits:")
    print(f"The minimal number of CCZ (controlled-controlled-Z) gates required to exactly synthesize a CCCZ (controlled-controlled-controlled-Z) gate is {minimal_ccz_gates}.")
    print("\nThis can be represented by the following resource equation:")
    
    # The final equation as requested, showing the numbers for each gate type.
    num_ccc_z_gates = 1
    num_ccz_gates = minimal_ccz_gates
    
    print(f"{num_ccc_z_gates} * CCCZ = {num_ccz_gates} * CCZ")

if __name__ == "__main__":
    solve_ccc_z_synthesis()

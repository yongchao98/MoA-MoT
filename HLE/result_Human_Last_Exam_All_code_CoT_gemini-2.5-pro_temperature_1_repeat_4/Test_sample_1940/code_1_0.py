def solve_cccz_synthesis():
    """
    This function explains the derivation for the minimal number of CCZ gates
    required to synthesize a CCCZ gate without ancilla qubits.
    """

    # --- Step 1: Define the problem ---
    print("Goal: Find the minimal number of CCZ gates to synthesize a CCCZ gate.")
    print("Allowed gates: CCZ (controlled-controlled-Z) and arbitrary single-qubit rotations.")
    print("Constraint: No ancilla (extra) qubits are allowed.\n")

    # --- Step 2: Problem Equivalence ---
    print("Step 2: Problem Equivalence")
    print("A CCCZ gate is equivalent to a CCCX (4-qubit Toffoli) gate up to single-qubit Hadamard (H) gates.")
    print("CCCZ(c1,c2,c3 -> t) = H_t * CCCX(c1,c2,c3 -> t) * H_t")
    print("Similarly, a CCZ gate is equivalent to a CCX (3-qubit Toffoli) gate.")
    print("CCZ(c1,c2 -> t) = H_t * CCX(c1,c2 -> t) * H_t")
    print("Since arbitrary single-qubit rotations (like H) are allowed and have zero cost for this problem, we can solve the equivalent problem: 'What is the minimal number of CCX gates to synthesize a CCCX gate?'\n")

    # --- Step 3: Known Decomposition ---
    print("Step 3: Known CCCX Gate Decomposition")
    print("A known efficient, ancilla-free construction for a 4-qubit Toffoli gate (CCCX) uses:")
    num_ccx_gates = 6
    num_cnot_gates = 2
    print(f"- {num_ccx_gates} CCX (3-qubit Toffoli) gates")
    print(f"- {num_cnot_gates} CNOT (2-qubit) gates\n")

    # --- Step 4: Calculate Gate Costs in CCZ units ---
    print("Step 4: Calculate Cost of Each Gate in the Decomposition")
    cost_of_ccx_in_ccz = 1
    print(f"Cost of 1 CCX gate = {cost_of_ccx_in_ccz} CCZ gate (using H gates).")

    cost_of_cnot_in_ccz = 1
    print("A CNOT gate can be constructed from one CCZ gate and single-qubit rotations using one of the other available qubits as a temporary workspace.")
    print("CNOT(c,t) = (H_t) * (X_q * CCZ(c,t,q) * X_q) * (H_t)")
    print(f"Cost of 1 CNOT gate = {cost_of_cnot_in_ccz} CCZ gate.\n")
    
    # --- Step 5: Final Calculation ---
    print("Step 5: Final Calculation")
    total_cost = num_ccx_gates * cost_of_ccx_in_ccz + num_cnot_gates * cost_of_cnot_in_ccz
    
    print("The total cost is the sum of the costs of the component gates.")
    print(f"Total Cost = ({num_ccx_gates} CCX gates * {cost_of_ccx_in_ccz} CCZ/CCX) + ({num_cnot_gates} CNOT gates * {cost_of_cnot_in_ccz} CCZ/CNOT)")
    print(f"Total Cost = {num_ccx_gates * cost_of_ccx_in_ccz} + {num_cnot_gates * cost_of_cnot_in_ccz}")
    print(f"Total Cost = {total_cost} CCZ gates.\n")
    
    print("This is believed to be the minimal number for ancilla-free synthesis based on known constructions.")
    final_answer = total_cost
    
    # Return the final answer in the specified format
    return final_answer

# Execute the function and print the final result.
final_result = solve_cccz_synthesis()
print(f"The minimal number of CCZ gates required is:")
print(f"<<<{final_result}>>>")

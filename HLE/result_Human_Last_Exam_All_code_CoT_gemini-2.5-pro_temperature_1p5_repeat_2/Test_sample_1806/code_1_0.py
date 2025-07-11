import numpy as np

def solve_quantum_stabilizer_problem():
    """
    Checks if a 4-qubit code is a stabilizer code for a given set of stabilizers.
    The code prints its reasoning and the result of the checks.
    """

    stabilizers = {
        "Z1 Z2": [('Z', 0), ('Z', 1)],
        "Z2 Z3": [('Z', 1), ('Z', 2)],
        "Z3 Z4": [('Z', 2), ('Z', 3)]
    }
    logical_states = {
        "|0_L>": 0,
        "|1_L>": 1
    }
    num_qubits = 4
    all_stabilized = True

    print("To be a stabilizer code, every logical basis state must be an eigenstate with eigenvalue +1 for every stabilizer generator.")
    print("We check this condition for S|psi_L> = eigenvalue * |psi_L>.\n")

    for s_name, s_ops in stabilizers.items():
        print(f"--- Checking Stabilizer {s_name} ---")
        for l_name, l_bit in logical_states.items():
            
            full_op_str = ['I'] * num_qubits
            for op_type, op_idx in s_ops:
                full_op_str[op_idx] = op_type
                
            term_eigenvalues = []
            for i in range(num_qubits):
                op = full_op_str[i]
                if op == 'Z':
                    # Z|0> = +1|0>, Z|1> = -1|1>
                    eigenvalue = 1 if l_bit == 0 else -1
                    term_eigenvalues.append(eigenvalue)
                else: # op == 'I'
                    term_eigenvalues.append(1)

            total_eigenvalue = np.prod(term_eigenvalues)

            # Build the string representation of the equation
            equation_lhs = ""
            for op_char, bit_val in zip(full_op_str, [l_bit]*num_qubits):
                equation_lhs += f"({op_char}|{bit_val}>)"
                
            eigenvalue_calc_str = ' * '.join(f"({v})" for v in term_eigenvalues)

            print(f"Applying {s_name} to {l_name}:")
            print(f"  {equation_lhs} = ({eigenvalue_calc_str}) |{str(l_bit)*num_qubits}>")
            print(f"  = {int(total_eigenvalue)} |{str(l_bit)*num_qubits}>")
            
            if total_eigenvalue != 1:
                all_stabilized = False
                print("  Result: NOT STABILIZED (Eigenvalue is not +1)\n")
            else:
                print("  Result: STABILIZED (Eigenvalue is +1)\n")

    print("--- Final Conclusion ---")
    if all_stabilized:
        print("Yes, the code can be considered a stabilizer code with the given stabilizers.")
        print("This is because all logical basis states are eigenvectors with eigenvalue +1 for all stabilizer generators.")
    else:
        print("No, the code cannot be considered a stabilizer code with the given stabilizers.")
        print("This is because at least one logical basis state is not an eigenvector with eigenvalue +1 for a stabilizer generator.")

solve_quantum_stabilizer_problem()
<<<Yes>>>
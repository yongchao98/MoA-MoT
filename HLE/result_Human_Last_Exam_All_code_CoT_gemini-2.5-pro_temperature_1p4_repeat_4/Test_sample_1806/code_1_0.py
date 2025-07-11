import numpy as np

def check_stabilizer_code():
    """
    Checks if the 4-qubit code defined by |0_L> = |0000> and |1_L> = |1111>
    is a stabilizer code with the given stabilizers.
    """
    
    # Define the logical states using their bit representations
    logical_states = {
        "|0_L>": "0000",
        "|1_L>": "1111"
    }

    # Define the stabilizers by the qubit indices they act on (0-indexed)
    stabilizers = {
        "S1 = Z1*Z2": [0, 1],
        "S2 = Z2*Z3": [1, 2],
        "S3 = Z3*Z4": [2, 3]
    }
    
    print("To check if the code is a stabilizer code with the given operators,")
    print("we verify if the logical basis states are eigenvectors with eigenvalue +1 for all stabilizers.\n")

    all_states_stabilized = True

    # Iterate through each logical state
    for state_name, state_bits in logical_states.items():
        print(f"--- Checking state {state_name} = |{state_bits}> ---")
        is_current_state_stabilized = True
        
        # Iterate through each stabilizer
        for stab_name, indices in stabilizers.items():
            qubit_idx1 = indices[0]
            qubit_idx2 = indices[1]
            
            # Get the bit values at the positions the stabilizer acts on
            bit1 = int(state_bits[qubit_idx1])
            bit2 = int(state_bits[qubit_idx2])
            
            # The eigenvalue for Z_i * Z_j on |...b_i...b_j...> is (-1)^(b_i + b_j)
            eigenvalue = (-1)**(bit1 + bit2)
            
            # Print the calculation steps
            print(f"Applying {stab_name}:")
            print(f"  The eigenvalue is calculated as (-1)^(bit{qubit_idx1+1} + bit{qubit_idx2+1}).")
            print(f"  For state {state_name}, this is (-1)^({bit1} + {bit2}).")
            
            sum_of_bits = bit1 + bit2
            print(f"  Equation: (-1)^({sum_of_bits}) = {eigenvalue}")
            
            if eigenvalue == 1:
                print("  Result: The eigenvalue is +1. The state is stabilized.\n")
            else:
                print(f"  Result: The eigenvalue is {eigenvalue}. The state is NOT stabilized.\n")
                is_current_state_stabilized = False
        
        if not is_current_state_stabilized:
            all_states_stabilized = False

    # Print the final conclusion
    print("--- Final Conclusion ---")
    if all_states_stabilized:
        print("Both logical basis states, |0_L> and |1_L>, are stabilized by all given operators.")
        print("Therefore, the code can be considered a stabilizer code with stabilizers Z1*Z2, Z2*Z3, and Z3*Z4.")
        final_answer = "Yes"
    else:
        print("At least one logical basis state is not stabilized by all operators.")
        print("Therefore, the code cannot be considered a stabilizer code with these stabilizers.")
        final_answer = "No"
    
    return final_answer

if __name__ == '__main__':
    answer = check_stabilizer_code()
    print(f"<<<{answer}>>>")

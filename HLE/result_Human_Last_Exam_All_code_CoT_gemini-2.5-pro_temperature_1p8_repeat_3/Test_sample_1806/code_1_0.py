import numpy as np

def check_stabilizer_code():
    """
    Checks if a 4-qubit code can be described by a given set of stabilizers.

    The code is defined by |0_L> = |0000> and |1_L> = |1111>.
    The stabilizers to check are S1=Z1Z2, S2=Z2Z3, and S3=Z3Z4.
    """
    # Define stabilizers as lists of Pauli operators (I, X, Y, Z)
    stabilizers = {
        "S1": ['Z', 'Z', 'I', 'I'],
        "S2": ['I', 'Z', 'Z', 'I'],
        "S3": ['I', 'I', 'Z', 'Z']
    }

    # Define the logical basis states
    logical_states = {
        "|0_L>": "0000",
        "|1_L>": "1111"
    }

    print("To be a valid stabilizer code, every logical state must be a +1 eigenstate of every stabilizer.")
    print("Let's check this condition: S|psi> = +1 * |psi>.\n")
    print("We recall that Z|0> = +1|0> and Z|1> = -1|1>.\n")
    
    all_stabilized = True

    for stab_name, stab_ops in stabilizers.items():
        # Create a printable string for the stabilizer, e.g., "Z1Z2"
        stab_str_print = "".join([f"Z{i+1}" for i, op in enumerate(stab_ops) if op == 'Z'])
        
        print(f"--- Checking Stabilizer {stab_name} = {stab_str_print} ---")

        for state_name, state_bits in logical_states.items():
            state_str_full = f"|{state_bits}>"
            eigenvalue = 1
            eigenvalue_calc_parts = []

            # Calculate the total eigenvalue and the breakdown of the calculation
            for i in range(len(state_bits)):
                op = stab_ops[i]
                qubit = state_bits[i]
                
                if op == 'Z':
                    if qubit == '0':
                        eigenvalue_calc_parts.append("(+1)")
                    else: # qubit == '1'
                        eigenvalue *= -1
                        eigenvalue_calc_parts.append("(-1)")
            
            # Format the calculation string, e.g., "(-1) * (-1)"
            eigenvalue_calc_str = " * ".join(eigenvalue_calc_parts)

            # Print the full equation
            print(f"{stab_str_print} {state_name} = {stab_str_print} {state_str_full} = {eigenvalue_calc_str} * {state_str_full} = {eigenvalue:+d} * {state_name}")

            if eigenvalue != 1:
                all_stabilized = False
        print() # Add a newline for better readability

    # Print the final conclusion
    print("--- Conclusion ---")
    if all_stabilized:
        print("Yes. All logical basis states are +1 eigenstates of all the proposed stabilizers.")
        print("The stabilizers also commute with each other. Therefore, the code can be considered a stabilizer code with stabilizers Z1Z2, Z2Z3, and Z3Z4.")
    else:
        print("No. At least one logical basis state is not a +1 eigenstate of a stabilizer.")

# Run the check
check_stabilizer_code()
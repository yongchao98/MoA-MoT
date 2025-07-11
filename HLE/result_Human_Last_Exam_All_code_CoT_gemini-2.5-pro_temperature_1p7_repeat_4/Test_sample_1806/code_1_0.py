import numpy as np

def check_stabilizer_action(stabilizer_name, stabilizer_qubits, state_name, state_bits):
    """
    Checks if a state is stabilized by a Z-type stabilizer and prints the equation.
    
    Args:
        stabilizer_name (str): The name of the stabilizer operator (e.g., "Z1*Z2").
        stabilizer_qubits (list): List of 0-indexed qubit indices the Z's act on.
        state_name (str): The name of the logical state (e.g., "|0_L>").
        state_bits (list): The list of bits representing the state (e.g., [0,0,0,0]).
        
    Returns:
        bool: True if the eigenvalue is +1, False otherwise.
    """
    # Calculate the exponent for the eigenvalue (-1)^exponent
    exponent = sum(state_bits[i] for i in stabilizer_qubits)
    eigenvalue = (-1)**exponent
    
    state_str = f"|{''.join(map(str, state_bits))}>"
    
    # Print the equation representing the stabilizer action
    print(f"  {stabilizer_name} {state_str} = (-1)^({'+'.join(str(state_bits[i]) for i in stabilizer_qubits)}) {state_str} = {eigenvalue:+d} * {state_str}")
    
    return eigenvalue == 1

def main():
    """
    Main function to check the stabilizer code conditions.
    """
    # Define the logical states in bit representation
    logical_states = {
        "|0_L>": [0, 0, 0, 0],
        "|1_L>": [1, 1, 1, 1]
    }
    
    # Define the stabilizers by their name and the qubits they act on (0-indexed)
    stabilizers = {
        "Z1*Z2": [0, 1],
        "Z2*Z3": [1, 2],
        "Z3*Z4": [2, 3]
    }
    
    print("Checking if the logical states are stabilized by the given operators...\n")
    
    all_stabilized = True
    
    # Iterate through all combinations of stabilizers and logical states
    for stab_name, stab_qubits in stabilizers.items():
        print(f"--- Applying Stabilizer {stab_name} ---")
        for state_name, state_bits in logical_states.items():
            if not check_stabilizer_action(stab_name, stab_qubits, state_name, state_bits):
                all_stabilized = False
        print()

    # Final conclusion
    print("--- Conclusion ---")
    if all_stabilized:
        print("Yes, this can be considered a stabilizer code with the given stabilizers.")
        print("All logical basis states are +1 eigenstates of all stabilizer operators.")
    else:
        print("No, this cannot be considered a stabilizer code with the given stabilizers.")
        print("At least one logical basis state was not a +1 eigenstate of a stabilizer.")

if __name__ == "__main__":
    main()

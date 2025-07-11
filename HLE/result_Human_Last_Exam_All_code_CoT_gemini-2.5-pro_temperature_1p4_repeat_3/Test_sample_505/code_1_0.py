def solve_qubit_calculation():
    """
    Calculates the maximum number of logical qubits for a given surface code configuration.
    
    The number of logical qubits in a single patch of surface code is equal to the number of holes.
    For multiple, independent patches, the total number of logical qubits is the sum of the qubits
    from each patch.
    """
    
    # Define the parameters of the problem
    num_patches = 2
    num_holes_per_patch = 2
    
    # Calculate logical qubits per patch
    logical_qubits_per_patch = num_holes_per_patch
    
    # Calculate the total number of logical qubits
    total_logical_qubits = num_patches * logical_qubits_per_patch
    
    # Print the explanation and the final equation
    print("The total number of logical qubits is calculated as follows:")
    print(f"Total Qubits = (Number of Patches) * (Number of Holes per Patch)")
    print(f"Total Qubits = {num_patches} * {num_holes_per_patch} = {total_logical_qubits}")

solve_qubit_calculation()
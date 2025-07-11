def check_stabilizer(stabilizer_name, stabilizer_indices, state_ket):
    """
    Checks if a state is a +1 eigenstate of a given stabilizer and prints the calculation.
    
    Args:
        stabilizer_name (str): The name of the stabilizer, e.g., "S1=Z1*Z2".
        stabilizer_indices (list): A list of 0-based qubit indices the stabilizer acts on.
        state_ket (str): A string representing the computational basis state, e.g., "0000".

    Returns:
        bool: True if the state is a +1 eigenstate, False otherwise.
    """
    print(f"Checking {stabilizer_name} on |{state_ket}>:")
    
    total_eigenvalue = 1
    eigenvalue_parts = []
    
    # Format the operator name part of the equation, e.g., "Z1*Z2"
    op_name = "*".join([f"Z{i+1}" for i in stabilizer_indices])

    # Calculate the total eigenvalue and the individual eigenvalue components for the equation
    for i in stabilizer_indices:
        if state_ket[i] == '1':
            total_eigenvalue *= -1
            eigenvalue_parts.append("(-1)")
        else: # state_ket[i] == '0'
            eigenvalue_parts.append("(+1)")

    # Construct and print the full equation demonstrating the action
    equation = f"{op_name} |{state_ket}> = "
    equation += " * ".join(eigenvalue_parts)
    equation += f" |{state_ket}> = ({total_eigenvalue})|{state_ket}>"
    print(equation)
    
    is_stabilized = (total_eigenvalue == 1)
    if is_stabilized:
        print(f"Result: The state |{state_ket}> is a +1 eigenstate of {stabilizer_name}.\n")
    else:
        print(f"Result: The state |{state_ket}> is NOT a +1 eigenstate of {stabilizer_name}.\n")
        
    return is_stabilized

def solve():
    """
    Main function to solve the problem.
    """
    logical_0 = '0000'
    logical_1 = '1111'
    stabilizers = {
        "S1=Z1*Z2": [0, 1],
        "S2=Z2*Z3": [1, 2],
        "S3=Z3*Z4": [2, 3]
    }

    print("To be a stabilizer code, the logical basis states must be +1 eigenstates of all stabilizer generators.")
    print("The generators must also commute, which is true here as they are all products of Z operators.")
    print("\nWe now check if the logical states are stabilized.\n")
    
    all_stabilized = True
    
    print("-" * 60)
    print(f"Checking logical state |0_L> = |{logical_0}>")
    print("-" * 60)
    for name, indices in stabilizers.items():
        if not check_stabilizer(name, indices, logical_0):
            all_stabilized = False
    
    print("-" * 60)
    print(f"Checking logical state |1_L> = |{logical_1}>")
    print("-" * 60)
    for name, indices in stabilizers.items():
        if not check_stabilizer(name, indices, logical_1):
            all_stabilized = False
            
    print("-" * 60)
    print("Final Conclusion:")
    if all_stabilized:
        print("Yes. Both |0_L> and |1_L> are +1 eigenstates of all stabilizer generators.")
        print("Therefore, the code can be considered a stabilizer code with the specified stabilizers.")
    else:
        print("No. One or more logical basis states are not +1 eigenstates of all the generators.")

solve()
<<<Yes>>>
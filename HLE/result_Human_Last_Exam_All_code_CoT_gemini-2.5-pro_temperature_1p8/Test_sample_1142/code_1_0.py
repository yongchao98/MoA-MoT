def solve_2d_chemistry():
    """
    This function calculates the 2D electron configuration for Nickel (Ni)
    and determines its bonding properties to answer the user's question.
    """
    atomic_number = 28
    element_name = "Ni"

    # The aufbau order of subshells in 2D
    subshell_order = ['1s', '2s', '2p', '3s', '3p', '3d', '4s', '4p', '4d', '5s']
    
    # Calculate capacities: s=2, others=4
    capacities = {s: 2 if s.endswith('s') else 4 for s in subshell_order}

    electrons_left = atomic_number
    configuration = []
    
    # Build the electron configuration string
    for subshell in subshell_order:
        if electrons_left == 0:
            break
        
        capacity = capacities[subshell]
        electrons_in_subshell = min(electrons_left, capacity)
        configuration.append(f"{subshell}{electrons_in_subshell}")
        electrons_left -= electrons_in_subshell
    
    # Get the last (valence) subshell and its properties
    last_subshell = configuration[-1]
    valence_shell_name = last_subshell[:-1]
    valence_electrons = int(last_subshell[-1])
    valence_capacity = capacities[valence_shell_name]
    
    # The final equation is the full electron configuration
    final_equation = " ".join(configuration)
    print(f"The 2D electron configuration for {element_name} (Z={atomic_number}) is:")
    print(final_equation)

    # Determine if it's a noble gas
    if valence_electrons == valence_capacity:
        conclusion = "The valence subshell is full, making the atom a noble gas."
        crystal_choice = "G"
        isotropy_answer = "n/a"
    else:
        # This case is not reached for Ni, but included for completeness
        conclusion = "The valence subshell is not full."
        crystal_choice = "Undetermined"
        isotropy_answer = "Undetermined"
        
    print(f"\nConclusion: {conclusion}")
    print("Therefore, no crystal structure forms.")

    # Print the final answer in the specified format for the two questions
    print(f"\nFinal Answer: {crystal_choice} {isotropy_answer}")


solve_2d_chemistry()
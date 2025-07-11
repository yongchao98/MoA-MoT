import math

def solve_constitution_game():
    """
    Calculates the smallest number of states in a "strong coalition"
    capable of ratifying a constitutional amendment against any opposition.
    """

    num_states = 50

    # Step 1: Calculate the number of states needed to force a proposal via a national convention.
    # This requires a 2/3 majority of states.
    states_for_proposal = math.ceil((2/3) * num_states)
    
    # Step 2: Calculate the number of states needed for ratification.
    # This requires a 3/4 majority of states.
    states_for_ratification = math.ceil((3/4) * num_states)
    
    # The size of the "strong coalition" must be sufficient for the most demanding step.
    min_coalition_size = max(states_for_proposal, states_for_ratification)
    
    # Print the explanation of the calculation.
    print("Analyzing the Constitutional Amendment process for a coalition vs. opposition:")
    print(f"Total number of states considered (N): {num_states}")
    print("-" * 30)
    
    print("Requirement 1: Proposing the Amendment")
    print("The only guaranteed path against opposition is calling a national convention.")
    print(f"This requires 2/3 of the states: ceil(2/3 * {num_states}) = {int(states_for_proposal)}")
    
    print("\nRequirement 2: Ratifying the Amendment")
    print("This requires 3/4 of the states.")
    print(f"This requires 3/4 of the states: ceil(3/4 * {num_states}) = {int(states_for_ratification)}")
    
    print("-" * 30)
    print("To be a 'strong coalition', the group must satisfy the larger of these two requirements.")
    final_equation = f"Smallest coalition size = max({int(states_for_proposal)}, {int(states_for_ratification)})"
    
    print(f"Final Equation: {final_equation}")
    print(f"Result: {min_coalition_size}")

solve_constitution_game()

import math

def solve_constitution_game():
    """
    Calculates the smallest number of states in a "strong coalition" capable
    of amending the U.S. Constitution.
    """
    num_states = 50

    # Fractions required by the Constitution
    convention_proposal_frac = 2/3
    ratification_frac = 3/4

    print("Analyzing the path to amend the Constitution via a National Convention.")
    print(f"The total number of states is {num_states}.")
    print("-" * 50)

    # Step 1: Calculate states needed to call a convention
    states_needed_for_convention = math.ceil(num_states * convention_proposal_frac)
    print("Step 1: Calling a National Convention")
    print("The Constitution requires 2/3 of state legislatures to call for a convention.")
    print(f"Calculation: ceil(2/3 * {num_states}) = ceil({num_states * convention_proposal_frac:.2f}) = {states_needed_for_convention} states.")
    print("-" * 50)

    # Step 2: Calculate states needed for ratification
    states_needed_for_ratification = math.ceil(num_states * ratification_frac)
    print("Step 2: Ratifying the Amendment")
    print("The Constitution requires 3/4 of states to ratify an amendment.")
    print(f"Calculation: ceil(3/4 * {num_states}) = ceil({num_states * ratification_frac:.2f}) = {states_needed_for_ratification} states.")
    print("-" * 50)

    # Step 3: Determine the minimum size of the strong coalition
    min_coalition_size = max(states_needed_for_convention, states_needed_for_ratification)
    print("Step 3: Determining the minimum size of a 'strong coalition'")
    print("A strong coalition must be large enough to both call the convention and ensure ratification.")
    print("Therefore, the required size is the larger of the two numbers.")
    print(f"Final Equation: Smallest Coalition Size = max({states_needed_for_convention}, {states_needed_for_ratification})")
    print(f"Result: The smallest number of states that could form a strong coalition is {min_coalition_size}.")

solve_constitution_game()